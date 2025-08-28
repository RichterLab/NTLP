import warnings

import numpy as np
import torch
import torch.nn as nn

from .data import read_training_file
from .physics import dydt,\
                     get_parameter_ranges, \
                     normalize_droplet_parameters, \
                     scale_droplet_parameters, \
                     solve_ivp_float32_outputs

class SimpleNet( nn.Module ):
    """
    4-layer multi-layer perceptron (MLP) with ReLU activations.  This aims to
    balance parameter count vs computational efficiency so that inferencing with
    it is faster than Gauss-Newton iterative solvers.
    """

    def __init__( self ):
        super().__init__()

        #
        # NOTE: These sizes were chosen without any consideration other than creating
        #       a small network (wrt parameter count) and should have good computational
        #       efficiency (wrt memory alignment and cache lines).  No effort has been
        #       spent to improve upon the initial guess.
        #
        self.fc1 = nn.Linear( 7, 32 )
        self.fc2 = nn.Linear( 32, 32 )
        self.fc3 = nn.Linear( 32, 32 )
        self.fc4 = nn.Linear( 32, 2 )

    def forward( self, x ):
        x = torch.relu( self.fc1( x ) )
        x = torch.relu( self.fc2( x ) )
        x = torch.relu( self.fc3( x ) )
        x = self.fc4( x )

        return x

class ResidualNet( SimpleNet ):
    """
    4-layer multi-layer perceptron (MLP) with ReLU activations.  This aims to
    balance parameter count vs computational efficiency so that inferencing with
    it is faster than Gauss-Newton iterative solvers.

    This residual network learns the delta between the input particle size and
    temperature and the outputs, given the provided background conditions.
    """

    def __init__( self ):
        super().__init__()

    def forward( self, x ):
        # Add the input to the final result to force the model to learn the
        # delta instead of the approximation itself.
        out  = super().forward( x )
        out += x[..., 0:2]

        return out

def do_bdf( input_parameters, times, gap_indices=np.array( [] ), **kwargs ):
    """
    Evaluates droplet parameters using backward differentiation formula (BDF) at the
    supplied simulation times.

    Evaluates each droplet parameter independently so as to provide comparison of
    BDF's solution against the reference observation.  The ith input is used to
    generate the ith output.

    Takes 4 arguments:

      input_parameters    - NumPy array, sized number_time_steps x 6, containing the
                            input parameters for a single particle in order by time.
                            These are provided in their natural, physical ranges.
      times               - NumPy array, shaped number_time_steps x 1, containing
                            the time at each step.
      gap_indices         - Optional NumPy array, shaped number_gaps x 1, containing
                            gap indices specifying the locations in input_parameters
                            immediately after a gap.  This is unused and only exists
                            to provide a compatible interface with do_iterative_bdf().
      kwargs              - Optional keyword arguments to supply to solve_ivp_float32_outputs().

    Returns 1 value:

      output_parameters   - NumPy array, sized number_time_steps x 2, containing
                            the evaluated size and temperature of the particle
                            along its background parameters with BDF in natural
                            ranges.

                            NOTE: The final timestep's outputs are zero since
                                  they cannot be computed without an additional
                                  simulation time.  This padding ensures
                                  output_parameters has a compatible shape with
                                  input_parameters.

    """

    #
    # NOTE: We provide gap_indices as an explicit argument to avoid manipulating
    #       our kwargs and removing it if it is present.  This routine does not
    #       care about gaps, but we generate a warning about an unused solver
    #       parameter in solve_ivp() if it is part of the kwargs passed on.
    #

    integration_times = np.diff( times )

    # Our outputs are sized the same as our inputs though the last value is set
    # to zero since we don't have an integration time (or output to compare
    # against) for the last observation.  As a result we skip it.
    output_parameters        = np.empty( (input_parameters.shape[0], 2), dtype=np.float32 )
    output_parameters[-1, :] = 0

    # Evaluate all but the last observation since those are the only ones we have
    # integration times for.
    for time_index in range( input_parameters.shape[0] - 1 ):
        output_parameters[time_index, :] = solve_ivp_float32_outputs( dydt,
                                                                      [0, integration_times[time_index]],
                                                                      input_parameters[time_index, :2],
                                                                      method="BDF",
                                                                      t_eval=[integration_times[time_index]],
                                                                      args=(input_parameters[time_index, 2:],),
                                                                      **kwargs )
        #
        # NOTE: We don't abort if this evaluation failed.  We hope that this
        #       particular observation was problematic but future observations
        #       aren't.
        #

    return output_parameters

def do_inference( input_parameters, times, model, device, **kwargs ):
    """
    Estimates droplet parameters using a specified model.  Model evaluation is
    performed on the device specified.

    Takes 5 arguments:

      input_parameters - NumPy array, sized number_time_steps x 6, containing
                         the input parameters for either one particle with
                         number_time_steps-many time steps or parameters for
                         number_time_steps-many droplets each with a single
                         observations.
      times            - NumPy array, sized number_time_steps x 1, containing
                         the time at each step.
      model            - PyTorch model to use.
      device           - String specifying the device to evaluate with.
      kwargs           - Optional keyword arguments.  These are provided for
                         interface compatibility with the other inference
                         routines and are ignored.

    Returns 1 value:

      output_parameters - NumPy array, sized number_droplets x 2, containing the
                          estimated radius and temperature for each of the droplets
                          at the specified integraition times.  These are in their
                          natural physical ranges.

                            NOTE: The final timestep's outputs are zero since
                                  they cannot be computed without an additional
                                  simulation time.  This padding ensures
                                  output_parameters has a compatible shape with
                                  input_parameters.

    """

    eval_model = model.to( device )
    eval_model.eval()

    # Stack the normalized parameters (all except the last) along with the
    # integration times so we have a single array to supply to the model.
    normalized_inputs = np.hstack( (normalize_droplet_parameters( input_parameters[:-1, :] ),
                                    np.diff( times ).reshape( (-1, 1) )) ).astype( "float32" )

    # Allocate the outputs and zero the last output since direct inference does
    # not produce one due to the lack of an integration time for the final
    # input.
    normalized_outputs        = torch.empty( (input_parameters.shape[0], 2), dtype=torch.float32 ).to( device )
    normalized_outputs[-1, :] = 0

    # Disable gradient accumulation since they're unnecessary for inference.
    with torch.no_grad():
        normalized_outputs[:-1, :] = eval_model( torch.from_numpy( normalized_inputs ).to( device ) )

    return scale_droplet_parameters( normalized_outputs.cpu() )

def do_iterative_bdf( input_parameters, times, gap_indices=np.array( [] ), **kwargs ):
    """
    Iteratively evaluates a particle trajectory along given background inputs with
    backward differentiation formula (BDF).  The outputs for the ith time are
    used as the input for the (i+1)th time.

    Takes 4 arguments:

      input_parameters    - NumPy array, sized number_time_steps x 6, containing the
                            input parameters for a single particle in order by time.
                            These are provided in their natural, physical ranges.
      times               - NumPy array, shaped number_time_steps x 1, containing
                            the time at each step.
      gap_indices         - Optional NumPy array, shaped number_gaps x 1, containing
                            gap indices specifying the locations in input_parameters
                            immediately after a gap.  May be supplied as an empty
                            array.  If omitted, gap detection is not performed,
                            otherwise the radius and temperatures from input_parameters
                            is used after a gap instead of the previously output
                            radius and temperature.
      kwargs              - Optional keyword arguments to supply to solve_ivp_float32_outputs().

    Returns 1 value:

      output_parameters   - NumPy array, sized number_time_steps x 2, containing the
                            estimated trajectory of a particle integrated
                            along its background parameters with BDF
                            in natural ranges.

                            NOTE: The first timestep's outputs are set to the
                                  inputs since they cannot be computed without
                                  an additional, preceding input.  This padding
                                  ensures output_parameters has a compatible
                                  shape with input_parameters.

    """

    integration_times = np.diff( times )

    # Construct a sequence of gap markers, including a sentinel that we'll never
    # reach, indicating when to reset our inputs.
    current_gap_number   = 0
    adjusted_gap_indices = np.append( gap_indices,
                                      [input_parameters.shape[0] + 1],
                                      axis=0 )

    output_parameters       = np.zeros( (input_parameters.shape[0], 2), dtype=np.float32 )
    output_parameters[0, :] = input_parameters[0, :2]

    # Evaluate all but the first observation, as it's our input.
    for time_index in range( 1, input_parameters.shape[0] ):
        # Skip evaluations which effectively don't take a step forward in time.
        if integration_times[time_index - 1] < 1.0e-10:
            output_parameters[time_index, :] = output_parameters[time_index - 1, :]
            continue

        # Select the correct input to this BDF evaluation based on whether we've
        # passed a gap or not.
        if time_index == adjusted_gap_indices[current_gap_number]:
            # Start again from the inputs when we've passed a gap.
            source_time_index  = time_index
            radius_temperature = input_parameters[source_time_index, :2]

            current_gap_number += 1
        else:
            # Otherwise continue iterating from the previous output.
            source_time_index  = time_index - 1
            radius_temperature = output_parameters[source_time_index, :2]

        # Solve for this time's output using the previous evaluation's output
        # as our input.
        output_parameters[time_index, :] = solve_ivp_float32_outputs( dydt,
                                                                      [0, integration_times[source_time_index]],
                                                                      radius_temperature,
                                                                      method="BDF",
                                                                      t_eval=[integration_times[source_time_index]],
                                                                      args=(input_parameters[source_time_index, 2:],),
                                                                      **kwargs )

        # See if this evaluation failed to provide a solution.  If so, then we
        # can short circuit the remaining evaluations since we don't have a way
        # to continue.
        if np.any( np.isnan( output_parameters[time_index, :] ) ):
            output_parameters[time_index:, :] = np.nan
            break

    return output_parameters

def do_iterative_inference( input_parameters, times, model, device, gap_indices=np.array( [] ) ):
    """
    Iteratively Estimates a particle trajectory along given background inputs
    using a MLP model.  The outputs for the ith time are used as the input for
    the (i+1)th time.

    Takes 5 arguments:

      input_parameters    - Array, sized number_time_steps x 6, containing the
                            input parameters for a single particle in order by time.
                            These are provided in their natural, physical ranges.
      times               - Array, sized number_time_steps x 1, containing the
                            time at each step.
      model               - PyTorch model to use.
      device              - String specifying the device to evaluate with.
      gap_indices         - Optional NumPy array, shaped number_gaps x 1, containing
                            gap indices specifying the locations in input_parameters
                            immediately after a gap.  May be supplied as an empty
                            array.  If omitted, gap detection is not performed,
                            otherwise the radius and temperatures from input_parameters
                            is used after a gap instead of the previously output
                            radius and temperature.

    Returns 1 value:

      output_parameters   - NumPy array, sized number_time_steps x 2, containing
                            the estimated trajectory of a particle integrated
                            along its background parameters with the MLP in
                            natural ranges.

                            NOTE: The first timestep's outputs are set to the
                                  inputs since they cannot be computed without
                                  an additional, preceding input.  This padding
                                  ensures output_parameters has a compatible
                                  shape with input_parameters.

    """

    eval_model = model.to( device )
    eval_model.eval()

    # Construct a sequence of gap markers, including a sentinel that we'll never
    # reach, indicating when to reset our inputs.
    current_gap_number   = 0
    adjusted_gap_indices = np.append( gap_indices,
                                      [input_parameters.shape[0] + 1],
                                      axis=0 )

    # Calculate integration time between time steps.  Keep an extra zero at the
    # end so that the array is the same length as the input parameters, though
    # only the first N-1 values are ever referenced.
    integration_times      = np.zeros( times.shape[0] )
    integration_times[:-1] = np.diff( times )

    normalized_data = torch.from_numpy( np.hstack( [normalize_droplet_parameters( input_parameters ),
                                                    integration_times.reshape( -1, 1 )] )
                                          .astype( "float32" ) ).to( device )

    with torch.no_grad():
        for time_index in range( 1,  normalized_data.shape[0] ):
            # Skip evaluations which effectively don't take a step forward in time.
            if integration_times[time_index - 1] < 1.0e-10:
                normalized_data[time_index, :2] = normalized_data[time_index - 1, :2]
                continue

            # Select the correct input to this BDF evaluation based on whether
            # we've passed a gap or not.
            if time_index == adjusted_gap_indices[current_gap_number]:
                # Start again from the inputs when we've passed a gap.
                source_time_index = time_index

                current_gap_number += 1
            else:
                # Otherwise continue iterating from the previous output.
                source_time_index = time_index - 1

            normalized_data[time_index, :2] = eval_model( normalized_data[source_time_index, :] )

    return scale_droplet_parameters( normalized_data[:, :2].to( "cpu" ).numpy() )

def generate_fortran_module( output_path, model_name, model, parameter_ranges={} ):
    """
    Creates a Fortran 2003 module that allows use of the supplied model with a
    batch size of 1 during inference.

    NOTE: This currently expects the supplied model state to represent a 4-layer
          MLP with a specific set weights/biases names.

    Takes 4 arguments:

      output_path      - File to write to.  If this exists it is overwritten.
      model_name       - Name of the model to write in the generated module's comments
                         so as to identify where the weights came from.
      model            - PyTorch model object for the module to expose.
      parameter_ranges - Optional dictionary of parameter ranges.  If omitted, defaults
                         to an empty dictionary and uses the current parameter ranges
                         returned by get_parameter_ranges().

    Returns nothing.

    """

    # Ensure we have a model we know how to serialize.
    if not isinstance( model, SimpleNet ):
        raise ValueError( "Provided model is type '{:s}' and not a subclass of SimpleNet!".format(
            str( type( model ) ) ) )

    #
    # NOTE: This is hard coded against the SimpleNet class' implementation.  This
    #       strange coupling could be better (less strange?) by moving this routine
    #       into SimpleNet.
    #
    expected_weights = ["fc1.weight",
                        "fc2.weight",
                        "fc3.weight",
                        "fc4.weight"]
    expected_biases  = ["fc1.bias",
                        "fc2.bias",
                        "fc3.bias",
                        "fc4.bias"]

    # Keep track of the parameter ranges for this model so we can report them at
    # run-time.
    if len( parameter_ranges ) == 0:
        model_parameter_ranges = get_parameter_ranges()
    else:
        model_parameter_ranges = parameter_ranges

    def _write_module_prolog( model_name, model_state, output_fp ):
        """
        Writes the module's prolog to the supplied file handle.

        Takes 3 arguments:

          model_name  - Name of the model to write in the generated module's comments
                        so as to identify where the weights came from.
          model_state - PyTorch model state dictionary for the model exposed by the module.
          output_fp   - File handle to write to.

        Returns nothing.

        """

        import datetime

        expected_parameters = set( expected_weights + expected_biases )

        # Confirm we have the weights and biases that we expect.
        if expected_parameters != set( model_state.keys() ):
            raise ValueError( "Model provided does not match a 4-layer MLP!" )

        # Ensure that the weights are 2D matrices.
        #
        # NOTE: We could be more thorough to ensure that the output of one layer
        #       matches the input of the next layer but we're in a rush right now...
        #
        for weights_name in expected_weights:
            if len( model_state[weights_name].shape ) != 2:
                raise ValueError( "'{:s}' is not a rank-2 tensor!".format( weights_name ) )

        # Ensure that the biases are vectors.
        #
        # NOTE: We could be more thorough to ensure that the output of one layer
        #       matches the input of the next layer but we're in a rush right now...
        #
        for bias_name in expected_biases:
            if len( model_state[bias_name].shape ) != 1:
                raise ValueError( "'{:s}' is not a rank-1 tensor!".format( bias_name ) )

        # Get the current date so people have an idea of when/how this module was
        # created.
        right_now             = datetime.datetime.now()
        current_date_time_str = right_now.strftime( "%Y/%m/%d at %H:%M:%S" )

        # Build a short string that describes the model's origins that can be
        # reported at run-time.
        model_metadata = "{:s} ({:s})".format(
            model_name,
            current_date_time_str )

        preamble_str = """!
! NOTE: This module was generated from {model_name:s} on {creation_date:s}.
!       Do not modify this directly!  Update the module generator's template and
!       regenerate this file!
!

module droplet_model

    ! Module containing a model of droplet size and temperature using a 4-layer
    ! multi-layer perceptron (MLP) neural network.  Once initialized, a single
    ! droplet's parameters and timestep can be used to estimate its future size
    ! and temperature.  This approach avoids an iterative search using
    ! Gauss-Newton to estimate its future charueacteristics and is typically much
    ! faster (O(30x) observed for small parameter count MLPs).
    !
    ! The module must be initialized with a call to initialize_model() and then
    ! droplet estimation, estimate(), can be performed as many times as
    ! necessary.  Weights for the MLP are hardcoded so initialization does not
    ! require any external resources.  An example calling sequence is the
    ! following:
    !
    !     ! Somewhere during simulation startup.
    !     call initialize_model()
    !
    !     ...
    !
    !     ! Particle simulation code.
    !     real*4 :: input_parameters(NUMBER_INPUTS)
    !     real*4 :: output_parameters(NUMBER_OUTPUTS)
    !
    !     do particle_index = 1, number_particles
    !
    !         input_parameters = [radius, temperature, salt_solute, air_temperature, rh, rhoa, t_final]
    !         call estimate( input_parameters, output_parameters )
    !
    !         new_radius      = output_parameters(RADIUS_INDEX)
    !         new_temperature = output_parameters(TEMPERATURE_INDEX)
    !
    !     end do
    !
    ! Inputs to the estimation, as are outputs, are provided in physical units
    ! and are internally normalized before using the model allowing this to be a
    ! drop-in replacement for alternative approaches for estimating droplet
    ! characteristics.
    !
    ! The model was trained on droplet parameters sampled from the product space
    ! of each of the individual parameters.  While this simple approach allowed
    ! for rapid development of a prototype it does include combinations of
    ! parameters that may be physically impossible which may have been
    ! under-sampled during training.  As a result the estimations in these
    ! corners of the parameter space may be of lower accuracy than more sampled,
    ! better behaved regions of the space.  The initial parameter ranges were
    ! selected by David Richter during prototyping and aim to cover all
    ! simulations of interest.
    !
    ! Things to note:
    !
    !   1. This approach was developed for code that would operate on a single
    !      droplet at a time vs multiple droplets at once.  As a result the
    !      application of a MLP network is simply a series of matrix-vector
    !      multiplications, vector-vector additions, and an activation function
    !      applied to the result.  We implement this case rather than attempt to
    !      generalize as the performance gain is already non-trivial and it is
    !      unclear when the code (if ever) would support batch estimation.
    !
    !   2. We could statically initialize each layer's weights and biases at
    !      compile time, albeit at the expense of readability.  By separating
    !      the weights and biases' definitions from their values we retain the
    !      ability to understand the structure and internal components of the
    !      MLP while requiring a single subroutine call at startup.
    !
    !   3. This is a naive implementation of a 4-layer MLP (it was written in 15
    !      minutes!)  and does not implement all potential optimizations.  In
    !      particular, each intermediate layer has its own array and no attempt
    !      to reuse arrays or work out of one large array has been made, so as
    !      to minimize the MLP's run-time footprint.  As initially developed the
    !      MLP is very small (O(2500) parameters) and uses O(10 KiB) storage,
    !      though should this ever be used for much larger MLPs, or with larger
    !      batch sizes (see above) then some of the omitted optimizations should
    !      be revisited.
    !
    !   4. While the MLP was trained on the product space of the input
    !      parameters and has poor performance on certain (likely) physically
    !      impossible parameter combinations, it is unclear what need to be done
    !      to improve the situation.  Re-training on a subset of the input
    !      parameters might improve performance and could be done without
    !      changing this implementation (albeit with undefined behavior for any
    !      physically impossible input) but it remains to be seen how much
    !      accuracy would be gained from this.
    !
    !   5. This approach is *NOT* OpenMP safe and is inherently single-threaded!
    !      Should multiple estimations need to be performed in parallel, an
    !      additional interface is needed that provides the intermediate layers
    !      to the estimation subroutine.  All of the weights and biases are
    !      effectively read-only and can be shared across threads.  This is
    !      straight forward to implement though it was not necessary for the
    !      initial implementation.
    !
    !   6. Single precision floating point is used throughout the module as it
    !      has been deemed to produce an acceptable level of accuracy for the
    !      radius and temperature parameters.  The input and output arrays
    !      must be real*4.
    !
    ! NOTE: This file was automatically generated from the the MLP's training
    !       code!  Do not modify this unless you really know what you're doing!
    !

    ! Model metadata to help with debugging.
    !
    ! NOTE: This is restricted to 128 characters to avoid running afoul of modern
    !       Fortran's 132 character limit.  We start the actual string literal on
    !       its own line to avoid generating a compiler warning when we have a
    !       long string.
    !
    character(len=128), parameter :: model_metadata = &
 "{model_metadata:s}"

    ! Indices into the input and output vectors.  The input vector uses all of
    ! the indices while the output vector only uses the first two.
    integer, parameter :: RADIUS_INDEX          = 1
    integer, parameter :: TEMPERATURE_INDEX     = 2
    integer, parameter :: SALT_SOLUTE_INDEX     = 3
    integer, parameter :: AIR_TEMPERATURE_INDEX = 4
    integer, parameter :: RH_INDEX              = 5
    integer, parameter :: RHOA_INDEX            = 6
    integer, parameter :: TFINAL_INDEX          = 7

    integer, parameter :: NUMBER_INPUTS  = {number_inputs:d}
    integer, parameter :: NUMBER_OUTPUTS = {number_outputs:d}

    integer, parameter :: NUMBER_HIDDEN_LAYER1_NEURONS = {number_layer1_neurons:d}
    integer, parameter :: NUMBER_HIDDEN_LAYER2_NEURONS = {number_layer2_neurons:d}
    integer, parameter :: NUMBER_HIDDEN_LAYER3_NEURONS = {number_layer3_neurons:d}

    ! The input variables must be normalized into the range of [-1, 1] before
    ! feeding them into the model.  The following ranges, means, and widths are
    ! used to map into, and from, the range of [-1, 1].
    !
    ! NOTE: These *must* match the training data!  Do not change these without
    !       retraining the model!
    !
    real*4, parameter :: RADIUS_LOG_RANGE(2)      = [{radius_start:.1f}, {radius_end:.1f}]
    real*4, parameter :: TEMPERATURE_RANGE(2)     = [{temperature_start:.1f}, {temperature_end:.1f}]
    real*4, parameter :: SALT_SOLUTE_LOG_RANGE(2) = [{salt_solute_start:.2f}, {salt_solute_end:.2f}]
    real*4, parameter :: AIR_TEMPERATURE_RANGE(2) = [{air_temperature_start:.1f}, {air_temperature_end:.1f}]
    real*4, parameter :: RH_RANGE(2)              = [{rh_start:.2f}, {rh_end:.2f}]
    real*4, parameter :: RHOA_RANGE(2)            = [{rhoa_start:.2f}, {rhoa_end:.2f}]

    real*4, parameter :: RADIUS_LOG_MEAN      = SUM( RADIUS_LOG_RANGE ) / 2
    real*4, parameter :: TEMPERATURE_MEAN     = SUM( TEMPERATURE_RANGE ) / 2
    real*4, parameter :: SALT_SOLUTE_LOG_MEAN = SUM( SALT_SOLUTE_LOG_RANGE ) / 2
    real*4, parameter :: AIR_TEMPERATURE_MEAN = SUM( AIR_TEMPERATURE_RANGE ) / 2
    real*4, parameter :: RH_MEAN              = SUM( RH_RANGE ) / 2
    real*4, parameter :: RHOA_MEAN            = SUM( RHOA_RANGE ) / 2

    !
    ! NOTE: We have be careful about lines no longer than 132 characters so we
    !       eliminate whitespace in the expressions.
    !
    real*4, parameter :: RADIUS_LOG_WIDTH      = (RADIUS_LOG_RANGE(2)-RADIUS_LOG_RANGE(1))/2
    real*4, parameter :: TEMPERATURE_WIDTH     = (TEMPERATURE_RANGE(2)-TEMPERATURE_RANGE(1))/2
    real*4, parameter :: SALT_SOLUTE_LOG_WIDTH = (SALT_SOLUTE_LOG_RANGE(2)-SALT_SOLUTE_LOG_RANGE(1))/2
    real*4, parameter :: AIR_TEMPERATURE_WIDTH = (AIR_TEMPERATURE_RANGE(2)-AIR_TEMPERATURE_RANGE(1))/2
    real*4, parameter :: RH_WIDTH              = (RH_RANGE(2)-RH_RANGE(1))/2
    real*4, parameter :: RHOA_WIDTH            = (RHOA_RANGE(2)-RHOA_RANGE(1))/2

    ! Weights for each of the layers.
    real*4, dimension(NUMBER_INPUTS, NUMBER_HIDDEN_LAYER1_NEURONS)                :: layer1_weights
    real*4, dimension(NUMBER_HIDDEN_LAYER1_NEURONS, NUMBER_HIDDEN_LAYER2_NEURONS) :: layer2_weights
    real*4, dimension(NUMBER_HIDDEN_LAYER2_NEURONS, NUMBER_HIDDEN_LAYER3_NEURONS) :: layer3_weights
    real*4, dimension(NUMBER_HIDDEN_LAYER3_NEURONS, NUMBER_OUTPUTS)               :: output_weights

    ! Biases for each of the layers.
    real*4, dimension(NUMBER_HIDDEN_LAYER1_NEURONS) :: layer1_biases
    real*4, dimension(NUMBER_HIDDEN_LAYER2_NEURONS) :: layer2_biases
    real*4, dimension(NUMBER_HIDDEN_LAYER3_NEURONS) :: layer3_biases
    real*4, dimension(NUMBER_OUTPUTS)               :: output_biases

    ! Arrays to hold the intermediate results between layers.
    !
    ! NOTE: We don't do anything fancy like detecting when we could reuse an
    !       intermediate.
    !
    real*4, dimension(NUMBER_HIDDEN_LAYER1_NEURONS) :: layer1_intermediate
    real*4, dimension(NUMBER_HIDDEN_LAYER2_NEURONS) :: layer2_intermediate
    real*4, dimension(NUMBER_HIDDEN_LAYER3_NEURONS) :: layer3_intermediate

    contains
""".format(
    # XXX: Double check the hidden layer neuron's sizes!  This was not
    #      thoroughly tested since the target model had identical sizes
    #      throughout.
    model_name=model_name,
    model_metadata=model_metadata,
    creation_date=current_date_time_str,
    number_inputs=model_state["fc1.weight"].shape[1],
    number_outputs=model_state["fc4.weight"].shape[0],
    number_layer1_neurons=model_state["fc1.weight"].shape[0],
    number_layer2_neurons=model_state["fc2.weight"].shape[0],
    number_layer3_neurons=model_state["fc3.weight"].shape[0],
    radius_start=model_parameter_ranges["radius"][0],
    radius_end=model_parameter_ranges["radius"][1],
    temperature_start=model_parameter_ranges["temperature"][0],
    temperature_end=model_parameter_ranges["temperature"][1],
    salt_solute_start=model_parameter_ranges["salt_solute"][0],
    salt_solute_end=model_parameter_ranges["salt_solute"][1],
    air_temperature_start=model_parameter_ranges["air_temperature"][0],
    air_temperature_end=model_parameter_ranges["air_temperature"][1],
    rh_start=model_parameter_ranges["relative_humidity"][0],
    rh_end=model_parameter_ranges["relative_humidity"][1],
    rhoa_start=model_parameter_ranges["rhoa"][0],
    rhoa_end=model_parameter_ranges["rhoa"][1]
        )

        print( preamble_str, file=output_fp )

    def _write_model_initializaiton_routine( model_state, output_fp, indentation_str ):
        """
        Writes the model's initialization routine to the supplied file handle with
        a specific amount of indentation.

        Takes 3 arguments:

          model_state     - PyTorch model state dictionary for the model to initialize.
          output_fp       - File handle to write to.
          indentation_str - String prepended to each line of the subroutine written.

        Returns nothing.

        """

        def _write_model_biases( model_state, output_fp, indentation_str ):
            """
            Internal routine that generates the biases initialization expressions for
            all of the biases in the supplied model state.  This does not generate
            initialization for weights nor does it generate a fully functional subroutine
            (i.e. pre-amble/epilog).

            Takes 3 arguments:

              model_state     - PyTorch model state dictionary for the model to initialize.
              output_fp       - File handle to write to.
              indentation_str - String prepended to each line of the subroutine written.

            Returns nothing.
            """

            # We rename each of the biases to Fortran-compatible symbols that are
            # self-descriptive.
            original_layer_names = expected_biases
            new_layer_names      = ["layer1_biases",
                                    "layer2_biases",
                                    "layer3_biases",
                                    "output_biases"]

            for original_layer_name, new_layer_name in zip( original_layer_names, new_layer_names ):
                layer_parameters = model_state[original_layer_name]

                # Start the array assignment to this bias' variable.
                biases_definition_str = "{:s}{:s} = [".format(
                    indentation_str,
                    new_layer_name )

                # Create enough indentation so all of the bias values are aligned
                # at the first column after the open bracket.
                biases_indentation_str = " " * len( biases_definition_str )

                # Template to join successive bias values together with.  Note that
                # this is not applied to the last one.
                definition_join_str    = ", &\n{:s}".format( biases_indentation_str )

                # Build the list of indented bias values, each printed with 10 digits
                # after the decimal point, and aligned so that positive biases have a
                # leading space to match negative biases' alignment.
                #
                # NOTE: This assumes the bias values' magnitudes are such that 10 digits
                #       of precision is sufficient.
                #
                bias_values_str = definition_join_str.join(
                    map( lambda x: "{: .10f}".format( x ),
                         layer_parameters.cpu().numpy() ) )

                print( "{:s}{:s}]".format(
                    biases_definition_str,
                    bias_values_str ),
                      file=output_fp )

        def _write_model_weights( model_state, output_fp, indentation_str ):
            """
            Internal routine that generates the weight initialization expressions for
            all of the weights in the supplied model state.  This does not generate
            initialization for biases nor does it generate a fully functional subroutine
            (i.e. pre-amble/epilog).

            Takes 3 arguments:

              model_state     - PyTorch model state dictionary for the model to initialize.
              output_fp       - File handle to write to.
              indentation_str - String prepended to each line of the subroutine written.

            Returns nothing.
            """

            # We rename each of the weights to Fortran-compatible symbols that are
            # self-descriptive.
            original_layer_names = expected_weights
            new_layer_names      = ["layer1_weights",
                                    "layer2_weights",
                                    "layer3_weights",
                                    "output_weights"]

            for original_layer_name, new_layer_name in zip( original_layer_names, new_layer_names ):
                layer_parameters = model_state[original_layer_name]

                # Get the shape of the array while reversing the order from
                # row-major (C, Python) to column-major (Fortran)
                outer_size, inner_size = layer_parameters.shape

                # Start the array assignment to this weight's variable.
                weights_definition_str = "{:s}{:s} = reshape( [".format(
                    indentation_str,
                    new_layer_name )

                # Create enough indentation so all of the weight's values are aligned
                # at the first column after the open bracket.
                weights_indentation_str = " " * len( weights_definition_str )

                # Template to join successive weight's values together with.  Note that
                # this is not applied to the last one.
                definition_join_str    = ", &\n{:s}".format( weights_indentation_str )

                # Build the list of indented weight values, each printed with 10 digits
                # after the decimal point, and aligned so that positive weights have a
                # leading space to match negative weights' alignment.  We ignore the
                # weight's shape since we'll reshape a 1D vector at run-time of the
                # compiled Fortran executable.
                #
                # NOTE: This assumes the weights values' magnitudes are such that 10 digits
                #       of precision is sufficient.
                #
                weights_values_str = definition_join_str.join(
                    map( lambda x: "{: .10f}".format( x ),
                         layer_parameters.cpu().numpy().ravel() ) )

                # Create the dimensions array supplied to reshape(), as well as the closing
                # parenthsis.
                weights_reshape_dimensions_str = "[{:d}, {:d}] )".format(
                    inner_size,
                    outer_size )

                print( "{:s}{:s}], &\n{:s}{:s}".format(
                    weights_definition_str,
                    weights_values_str,
                    weights_indentation_str,
                    weights_reshape_dimensions_str ),
                      file=output_fp )

        def _write_model_metadata( model_parameter_ranges, output_fp, indentation_str ):
            """
            Internal routine that generates code for the root MPI rank to report
            information about the MLP model.  This prints the Fortran
            droplet_model's model_metadata variable as well as pretty prints the
            model's parameter ranges to standard output so it is available in a
            simulation's log.  These aim to improve debugging issues with MLP
            inference.

            Takes 3 arguments:

              model_parameter_ranges - Dictionary of parameter ranges the model was trained on.
                                       See get_parameter_ranges() for details.
              output_fp              - File handle to write to.
              indentation_str        - String prepended to each line of the subroutine written.

            Returns nothing.
            """

            # Report a pretty printed table of parameter ranges along with their units.
            #
            # NOTE: The lack of indentation is intentional as these lines are long and we attempt
            #       to keep them under the 132 character limit by skipping indentation inside of
            #       the rank check conditional.
            #
            print_ranges_str = \
"""
! Write the model's metadata and pretty print a table of parameter
! ranges and their units to the log.
if( myid == 0 ) then
write(*, "(A,A,/)" ) "MLP model metadata: ", model_metadata

write(*, "(/,A,/,/,A,F5.2,A,F5.2,A,/,A,F6.2,A,F6.2,A,/,A,F6.2,A,F6.2,A)" ) &
                     "MLP was trained with the following parameter ranges:", &
                     "    Log10(Radius):      [ ", RADIUS_LOG_RANGE(1), ",  ", RADIUS_LOG_RANGE(2), "] m", &
                     "    Temperature:        [", TEMPERATURE_RANGE(1), ", ", TEMPERATURE_RANGE(2), "] K", &
                     "    log10(Salt solute): [", SALT_SOLUTE_LOG_RANGE(1), ", ", SALT_SOLUTE_LOG_RANGE(2), "] kg"
write(*, "(A,F6.2,A,F6.2,A,/,A,F6.2,A,F6.2,A,/,A,F4.2,A,F4.2,A,/)" ) &
                     "    Air temperatures:   [", AIR_TEMPERATURE_RANGE(1), ", ", AIR_TEMPERATURE_RANGE(2), "] K", &
                     "    Relative humidity:  [", RH_RANGE(1) * 100.0, ", ", RH_RANGE(2) * 100.0, "] %", &
                     "    Air density:        [  ", RHOA_RANGE(1), ",   ", RHOA_RANGE(2), "] kg/m^3"
end if
"""

            for print_ranges_line in print_ranges_str.splitlines():
                print( "{:s}{:s}".format(
                    indentation_str if len( print_ranges_line ) > 0 else "",
                    print_ranges_line ),
                       file=output_fp )

        # Our initialization routine is inside of a "contains" block so we're
        # indented one level deeper than the caller.
        inner_indentation_str = "    {:s}".format( indentation_str )

        # Prolog and epilog of our initialization routine, along with extra
        # newlines so we survive the application of .splitlines() below.
        #
        # NOTE: We need the current MPI rank from the pars module so we can
        #       properly report our state from a single rank.
        #
        initialization_prolog_str = ("subroutine initialize_model()\n" +
                                     "\n" +
                                     "    use pars, only:    myid\n" +
                                     "\n"
                                     "    implicit none\n" +
                                     "\n")
        initialization_epilog_str = ("\n" +
                                     "end subroutine initialize_model")

        # Write the subroutine's prolog.
        for initialization_line in initialization_prolog_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( initialization_line ) > 0 else "",
                initialization_line ),
                file=output_fp )

        # Write the weights and biases separated by an empty line.
        _write_model_weights( model_state, output_fp, inner_indentation_str )
        print( "", file=output_fp )
        _write_model_biases( model_state, output_fp, inner_indentation_str )

        # Report the parameter ranges at run-time so users have an idea of where
        # the model should perform well.
        print( "", file=output_fp )
        _write_model_metadata( model_parameter_ranges, output_fp, inner_indentation_str )

        # Write the subroutine's epilog.
        for initialization_line in initialization_epilog_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( initialization_line ) > 0 else "",
                initialization_line ),
                file=output_fp )

    def _write_model_inference_routine( model_state, output_fp, indentation_str ):
        """
        Writes the model's inference routine to the supplied file handle with
        a specific amount of indentation.

        Takes 3 arguments:

          model_state     - Unused argument kept for a consistent signature.
          output_fp       - File handle to write to.
          indentation_str - String prepended to each line of the subroutine written.

        Returns nothing.

        """

        residual_inference_str = """
    ! Add the model's input to its output since it's learned to compute a delta rather
    ! than the final answer.
    output(RADIUS_INDEX)      = output(RADIUS_INDEX)      + normalized_input(RADIUS_INDEX)
    output(TEMPERATURE_INDEX) = output(TEMPERATURE_INDEX) + normalized_input(TEMPERATURE_INDEX)
"""

        # Only include the residual if the model supports it.
        if not isinstance( model, ResidualNet ):
            residual_inference_str = ""

        inference_str = """
! Estimates a droplet's future size and temperature based on it's current size
! and temperature, as well as key parameters from the environment its in.
! Inputs and outputs are provided/returned in un-normalized physical scales.
subroutine estimate( input, output )

    implicit none

    !
    ! NOTE: This is hardcoded to assume a batch size of 1 so we can have
    !       a 1D array vs dealing with a 2D array with a singleton
    !       dimension.
    !
    real*4, dimension(NUMBER_INPUTS), intent(in)   :: input
    real*4, dimension(NUMBER_OUTPUTS), intent(out) :: output

    real*4, dimension(NUMBER_INPUTS)               :: normalized_input
    integer                                        :: output_index

    ! Normalize the non-temporal inputs so they're in the range [-1, 1].
    normalized_input(RADIUS_INDEX)          = (log10( input(RADIUS_INDEX) ) - RADIUS_LOG_MEAN) / RADIUS_LOG_WIDTH
    normalized_input(TEMPERATURE_INDEX)     = (input(TEMPERATURE_INDEX) - TEMPERATURE_MEAN) / TEMPERATURE_WIDTH
    normalized_input(SALT_SOLUTE_INDEX)     = (log10( input(SALT_SOLUTE_INDEX) ) - SALT_SOLUTE_LOG_MEAN) / SALT_SOLUTE_LOG_WIDTH
    normalized_input(AIR_TEMPERATURE_INDEX) = (input(AIR_TEMPERATURE_INDEX) - AIR_TEMPERATURE_MEAN) / AIR_TEMPERATURE_WIDTH
    normalized_input(RH_INDEX)              = (input(RH_INDEX) - RH_MEAN) / RH_WIDTH
    normalized_input(RHOA_INDEX)            = (input(RHOA_INDEX) - RHOA_MEAN) / RHOA_WIDTH

    ! Our integration time remains as is.
    normalized_input(TFINAL_INDEX)          = input(TFINAL_INDEX)

    ! Compute x_1 = ReLU( W_1*I + b_1 ).
    do output_index = 1, NUMBER_HIDDEN_LAYER1_NEURONS
        layer1_intermediate(output_index) = &
             max( sum( normalized_input(:) * layer1_weights(:, output_index) ) + layer1_biases(output_index), 0.0 )
    end do

    ! Compute x_2 = ReLU( W_2*x_1 + b_2 ).
    do output_index = 1, NUMBER_HIDDEN_LAYER2_NEURONS
        layer2_intermediate(output_index) = &
             max( sum( layer1_intermediate(:) * layer2_weights(:, output_index) ) + layer2_biases(output_index), 0.0 )
    end do

    ! Compute x_3 = ReLU( W_3*x_2 + b_3 ).
    do output_index = 1, NUMBER_HIDDEN_LAYER3_NEURONS
        layer3_intermediate(output_index) = &
             max( sum( layer2_intermediate(:) * layer3_weights(:, output_index) ) + layer3_biases(output_index), 0.0 )
    end do

    ! Compute O = W_4*x_3 + b_4.
    do output_index = 1, NUMBER_OUTPUTS
        output(output_index) = sum( layer3_intermediate(:) * output_weights(:, output_index) ) + output_biases(output_index)
    end do
{residual_inference:s}
    ! Scale the outputs to the expected ranges.
    output(RADIUS_INDEX)      = 10.0**(output(RADIUS_INDEX) * RADIUS_LOG_WIDTH + RADIUS_LOG_MEAN)
    output(TEMPERATURE_INDEX) = output(TEMPERATURE_INDEX) * TEMPERATURE_WIDTH + TEMPERATURE_MEAN

end subroutine estimate

""".format(
    residual_inference=residual_inference_str
)

        for inference_line in inference_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( inference_line ) > 0 else "",
                inference_line ),
                   file=output_fp )

    def _write_module_epilog( model_state, output_fp ):
        """
        Writes the module's epilog to the supplied file handle.

        Takes 2 arguments:

          model_state     - Unused argument kept for a consistent signature.
          output_fp       - File handle to write to.

        Returns nothing.

        """

        epilog_str = "end module droplet_model"

        print( epilog_str, file=output_fp )

    # Module subroutines, inside of a "contains" block, are indented twice.
    indentation_str = " " * 8

    # Convenience variable for the internal methods that take the model weights
    # instead of the model itself.
    model_state = model.state_dict()

    with open( output_path, "w" ) as output_fp:
        # Write out the beginning of the module.  This includes the definitions for
        # the layers' weights and biases but not the method that initializes them.
        _write_module_prolog( model_name, model_state, output_fp )

        # Write out both of the subroutines that comprise this model.
        _write_model_initializaiton_routine( model_state, output_fp, indentation_str )
        _write_model_inference_routine( model_state, output_fp, indentation_str )

        # Write out the end of the module.
        _write_module_epilog( model_state, output_fp )

def load_model_checkpoint( checkpoint_path, model, optimizer=None ):
    """
    Loads a checkpoint from disk and updates the provided model and optimizer
    objects with the state contained.  Returns the parameter ranges found in the
    checkpoint as well as the model's recorded training loss.

    Takes 3 arguments:

      checkpoint_path - Path to the file that checkpoint was loaded from.
      model           - Torch model whose weights will be set to the
                        checkpoint's contents.
      optimizer       - Optional Torch optimizer whose state will be set to the
                        checkpoint's contents.  If omitted, any optimizer state
                        in checkpoint_path is ignored.

    Returns 3 values:

      parameter_ranges - Dictionary containing the droplet parameter ranges
                         associated with model.  See get_parameter_range() for
                         details.
      loss_function    - Function handle for the loss function associated with
                         model.
      training_loss    - Sequence containing model's training loss across
                         batches.

    """

    def _load_model_checkpoint_v1( checkpoint_path, model, checkpoint ):
        """
        Loads a legacy, version 1 checkpoint.  This is simply a serialization of
        the model's state dictionary without optimizer state, the loss function,
        or training lost history.  The model provided is updated in-place.

        Takes 3 arguments:

          checkpoint_path - Path to the file that checkpoint was loaded from.
          model           - Torch model whose weights will be set to the
                            checkpoint's contents.
          checkpoint      - Dictionary containing the weights and biases for
                            model.  This must be suitable for loading via
                            .load_state_dict().

        Returns 1 value:

          parameter_ranges - Dictionary containing the current droplet parameter
                             ranges.  A warning is reported since it is not
                             guaranteed that these ranges match the model
                             itself.  See get_parameter_range() for details on
                             the dictionary's keys.

        """

        # Legacy checkpoints do not contain parameter ranges so warn the user
        # that what we're returning may not match the model that was loaded.
        warnings.warn( "Loaded a legacy checkpoint from '{:s}', returning current parameter ranges!".format(
            checkpoint_path ) )

        warnings.warn( "Version 1 checkpoints are incompatible with the current dydt()!  "
                       "Make sure to convert salt solute to salt mass before inferencing." )

        # The checkpoint is simply the model's weights and biases.
        model.load_state_dict( checkpoint )

        return get_parameter_ranges()

    def _load_model_checkpoint_v2( checkpoint_path, model, optimizer, checkpoint ):
        """
        Loads a version 2 checkpoint.  These contain model weights, the droplet
        parameter ranges the model was trained on, and the optimizer's state.
        The model and optimizer are updated in-place.  The model and parameter
        ranges use salt mass.

        Takes 4 arguments:

          checkpoint_path - Path to the file that checkpoint was loaded from.
          model           - Torch model whose weights will be set to the
                            checkpoint's contents.
          optimizer       - Optional Torch optimizer whose state will be set to
                            the checkpoint's contents.  If omitted, the loaded
                            optimizer state is discarded.
          checkpoint      - Dictionary containing the following keys:
                            "droplet_parameters", "loss_function",
                            "model_weights", "optimizer_state".

        Returns 3 values:

          parameter_ranges - Dictionary containing the droplet parameter ranges
                             associated with model.  See get_parameter_range()
                             for details.
          loss_function    - Function handle for the loss function associated
                             with model.
          training_loss    - Sequence containing model's training loss across
                             batches.

        """

        model.load_state_dict( checkpoint["model_weights"] )
        if optimizer is not None:
            optimizer.load_state_dict( checkpoint["optimizer_state"] )

        parameter_ranges = checkpoint["droplet_parameter_ranges"]
        loss_function    = checkpoint["loss_function"]
        training_loss    = checkpoint["training_loss"]

        # Sanity check that the parameter ranges are as expected.
        if "salt_mass" not in parameter_ranges:
            raise ValueError( "'{:s}' is not a valid version 2 checkpoint.  "
                              "Its parameter ranges does not contain 'salt mass'!".format(
                                  checkpoint_path ) )

        warnings.warn( "Version 2 checkpoints are incompatible with the current dydt()!  "
                       "Make sure to convert salt solute to salt mass before inferencing." )

        return parameter_ranges, loss_function, training_loss

    def _load_model_checkpoint_v3( checkpoint_path, model, optimizer, checkpoint ):
        """
        Loads a version 3 checkpoint.  These contain model weights, the droplet
        parameter ranges the model was trained on, and the optimizer's state.
        The model and optimizer are updated in-place.  The model and parameter
        ranges use salt solute instead of salt mass.

        Takes 4 arguments:

          checkpoint_path - Path to the file that checkpoint was loaded from.
          model           - Torch model whose weights will be set to the
                            checkpoint's contents.
          optimizer       - Optional Torch optimizer whose state will be set to
                            the checkpoint's contents.  If omitted, the loaded
                            optimizer state is discarded.
          checkpoint      - Dictionary containing the following keys:
                            "droplet_parameters", "loss_function",
                            "model_weights", "optimizer_state".

        Returns 3 values:

          parameter_ranges - Dictionary containing the droplet parameter ranges
                             associated with model.  See get_parameter_range()
                             for details.
          loss_function    - Function handle for the loss function associated
                             with model.
          training_loss    - Sequence containing model's training loss across
                             batches.

        """

        model.load_state_dict( checkpoint["model_weights"] )
        if optimizer is not None:
            optimizer.load_state_dict( checkpoint["optimizer_state"] )

        parameter_ranges = checkpoint["droplet_parameter_ranges"]
        loss_function    = checkpoint["loss_function"]
        training_loss    = checkpoint["training_loss"]

        # Sanity check that the parameter ranges are as expected.
        if "salt_solute" not in parameter_ranges:
            print( parameter_ranges )
            raise ValueError( "'{:s}' is not a valid v3 checkpoint.  "
                              "Its parameter ranges does not contain 'salt solute'!".format(
                                  checkpoint_path ) )

        return parameter_ranges, loss_function, training_loss

    # Load the checkpoint's Tensors as a dictionary and figure out which version
    # it is.
    checkpoint         = torch.load( checkpoint_path, weights_only=False, map_location=torch.device('cpu') )
    checkpoint_version = checkpoint.get( "checkpoint_version", 1 )

    # Load the checkpoint based on its reported version.  Complain if we don't
    # know how to support it.
    if checkpoint_version == 1:
        # Legacy checkpoints did not track the loss function or training
        # loss.
        loss_function = None
        training_loss = []

        parameter_ranges = _load_model_checkpoint_v1( checkpoint_path,
                                                      model,
                                                      checkpoint )
    elif checkpoint_version == 2:
        (parameter_ranges,
         loss_function,
         training_loss) = _load_model_checkpoint_v2( checkpoint_path,
                                                     model,
                                                     optimizer,
                                                     checkpoint )
    elif checkpoint_version == 3:
        (parameter_ranges,
         loss_function,
         training_loss) = _load_model_checkpoint_v3( checkpoint_path,
                                                     model,
                                                     optimizer,
                                                     checkpoint )
    else:
        raise RuntimeError( "Unknown checkpoint version ({:d}) in '{:s}'!".format(
            checkpoint_version,
            checkpoint_path ) )

    return parameter_ranges, loss_function, training_loss

def ode_residual( inputs, outputs, model ):

    parameter_ranges = get_parameter_ranges()

    drdt = torch.autograd.grad( outputs[:, 0],
                                inputs,
                                grad_outputs=torch.ones_like( outputs[:, 0] ),
                                create_graph=True )[0]
    dTdt = torch.autograd.grad( outputs[:, 1],
                                inputs,
                                grad_outputs=torch.ones_like( outputs[:, 1] ),
                                create_graph=True )[0]

    drdt *= np.diff( parameter_ranges["radius"] ).astype( float ) * (10**((outputs[:, 0] * np.diff( parameter_ranges["radius"] ).astype( float ) / 2) + np.mean( parameter_ranges["radius"] ).astype( float ))) * 0.5 * np.log( 10 )
    dTdt *= np.diff( parameter_ranges["temperature"] ).astype( float ) * 0.5

    return [drdt, dTdt]

def save_model_checkpoint( checkpoint_prefix, checkpoint_number, model, optimizer, loss_function, training_loss, parameter_ranges={} ):
    """
    Serializes the provided model and optimizer state disk so that it may
    be loaded for additional training and/or inferencing.  Existing checkpoints
    are overwritten.

    Takes 7 arguments:

      checkpoint_prefix - Prefix of the checkpoint path.
      checkpoint_number - Number of the checkpoint to save.  When non-negative
                          is combined with checkpoint_prefix to create checkpoint_path.
      model             - Torch model object whose state should be saved.
      optimizer         - Torch optimizer object whose state should be saved.
      loss_function     - Function that computes the loss for model.
      training_loss     - Sequence containing the training loss history.
      parameter_ranges  - Optional dictionary containing the droplet parameter
                          ranges used to train model.  If omitted, defaults to
                          an empty dictionary and the current parameter ranges
                          are stored in the checkpoint.

    Returns 1 value:

      checkpoint_path - Path to the checkpoint written, with a ".pt" extension.
                        If checkpoint_number is negative then checkpoint_prefix
                        simply has an extension appended, otherwise
                        checkpoint_number is added before the extension.

    """

    #
    # NOTE: We only support saving models in the most recent version.  There
    #       isn't a need to be fully flexible during development and one-off
    #       cases where we need to write a specific version can be handled
    #       as needed.
    #
    CHECKPOINT_VERSION = 3

    # Build the checkpoint path.
    if checkpoint_number >= 0:
        checkpoint_path = "{:s}_{:d}.pt".format(
            checkpoint_prefix,
            checkpoint_number )
    else:
        checkpoint_path = "{:s}.pt".format(
            checkpoint_prefix )

    current_parameter_ranges = get_parameter_ranges()

    # Update the parameter ranges with the caller's.
    if len( parameter_ranges ) > 0:
        current_parameter_ranges.update( parameter_ranges )

    # Checkpoints are dictionaries comprised of the following key/values:
    #
    #  droplet_parameter_ranges - The physics module's current parameter range
    #                             dictionary.
    #  loss_function            - Function that computes model's training loss.
    #  model_weights            - Model weights and biases.
    #  optimizer_state          - Optional optimizer state so as to continue
    #                             training where the checkpoint left off.
    #  training_loss            - Sequence of training loss values.
    #
    checkpoint = {
        "checkpoint_version":       CHECKPOINT_VERSION,
        "droplet_parameter_ranges": current_parameter_ranges,
        "loss_function":            loss_function,
        "model_weights":            model.to( "cpu" ).state_dict(),
        "optimizer_state":          optimizer.state_dict(),
        "training_loss":            training_loss
    }

    # Serialize the checkpoint.
    torch.save( checkpoint, checkpoint_path )

    return checkpoint_path

def train_model( model, criterion, optimizer, device, number_epochs, training_file, validation_file=None, checkpoint_prefix=None, epoch_callback=None, lr_scale=0.5, batch_size=1024 ):
    """
    Trains the supplied model for one or more epochs using all of the droplet parameters
    in an on-disk training file.  The parameters are read into memory once and then
    randomly sampled each epoch.  Any weird parameters encountered in the training file
    are logged when requested.  Model performance may be evaluated when a validation
    file is provided and is done so at the end of each epoch.

    Takes 11 arguments:

      model             - PyTorch model to optimize.
      criterion         - PyTorch loss object to use during optimization.
      optimizer         - PyTorch optimizer associated with model.
      device            - Device string indicating where the optimization is
                          being performed.
      number_epochs     - Number of epochs to train model for.  All training
                          data in training_file will be seen by the model
                          this many times.
      training_file     - Path to the file containing training data created by
                          create_training_file() OR a tuple containing three NumPy
                          arrays: input parameters, output parameters, and
                          integration times (sized N x 6, N x 2, and N x 1,
                          respectively, where N are the number of training
                          samples).
      validation_file   - Optional path to the file containing validation data created by
                          create_training_file() OR a tuple containing three NumPy
                          arrays: input parameters, output parameters, and
                          integration times (sized N x 6, N x 2, and N x 1,
                          respectively, where N are the number of training
                          samples).  If omitted, defaults to None and no
                          validation evaluations are performed
      checkpoint_prefix - Optional path prefix where model checkpoints should
                          be stored.  If omitted, defaults to None and
                          checkpoints are not written.  Otherwise, the epoch
                          number is appended to construct the path where each
                          epoch's checkpoint is written.
      epoch_callback    - Optional function to be called after each training
                          epoch.  If omitted, defaults to None.  When provided
                          must take:

                            model:           Torch model object
                            epoch_number:    Training epoch that just completed
                            optimizer:       Torch optimizer object
                            training_loss:   Sequence of training loss for the
                                             last epoch
                            validation_loss: Scalar validation loss for the last
                                             epoch

      lr_scale          - Optional floating point value, in the range of (0, 1],
                          to scale the learning rate at the end of each epoch.
                          If omitted, defaults to 0.5.
      batch_size        - Optional batch size used during training.  If omitted,
                          defaults to 1024 samples.

    Returns 2 values:

      training_loss_history   - List of training losses, one per mini-batch, per
                                epoch, during the training process.
      validation_loss_history - List of validation losses, one per epoch, during
                                the training process.

    """

    # Number of batches to train on per training loss measurement.
    MINI_BATCH_SIZE = 1024

    # Number of validation samples per validation batch.
    VALIDATION_BATCH_SIZE = 10240

    # Setup our callback(s) so we simply call them at the right time during
    # training.
    epoch_callbacks = []
    if epoch_callback is not None:

        # Support being provided a single callback instead of a list of them.
        if callable( epoch_callback ):
            epoch_callbacks.append( epoch_callback )
        # Blow up if we didn't have a callback or a list of them.
        elif not isinstance( epoch_callback, list ):
            raise ValueError( "Epoch callback is neither callable nor a list ({:s})!".format(
                type( epoch_callback ) ) )
        elif not all( map( lambda x: callable( x ), epoch_callback ) ):
            raise ValueError( "One of the callbacks provided isn't callable!" )
        else:
            # Add all of the callbacks.
            epoch_callbacks.extend( epoch_callback )

    #
    # NOTE: This is inefficient and requires the entire training data set to reside in
    #       RAM.  We could read chunks of the file on demand but that would require
    #       a more sophisticated training loop that performs I/O in a separate thread
    #       while the optimization process executes.
    #
    #       TL;DR Find a big machine to train on.
    #
    if isinstance( training_file, tuple ):
        training_inputs, training_outputs, training_times = training_file
    else:
        training_inputs, training_outputs, training_times = read_training_file( training_file )

    if validation_file is not None:
        if isinstance( validation_file, tuple ):
            validation_inputs, validation_outputs, validation_times = validation_file
        else:
            validation_inputs, validation_outputs, validation_times = read_training_file( validation_file )

        # Normalize the validation data once and prepare it for inference.
        #
        # NOTE: We don't send it to the device since it may be too large
        #       to hold on to during training.
        #
        validation_inputs  = normalize_droplet_parameters( validation_inputs )
        validation_outputs = normalize_droplet_parameters( validation_outputs )

        validation_inputs  = np.hstack( (validation_inputs,
                                         validation_times.reshape( -1, 1 )) )

        validation_inputs  = torch.from_numpy( validation_inputs )
        validation_outputs = torch.from_numpy( validation_outputs )

    weights = np.reciprocal( training_times )
    weights = np.stack( (weights, weights), axis=-1 )
    weights = torch.from_numpy( weights ).to( device )

    #
    # NOTE: This ignores parameters if the last batch isn't complete.
    #
    number_batches = training_inputs.shape[0] // batch_size

    # Track each mini-batch's training loss for analysis.
    training_loss_history = []

    # Track each epoch's validation loss for analysis.
    validation_loss_history = []

    for epoch_index in range( number_epochs ):
        model.train()

        # Reset our training loss.
        training_loss = []
        running_loss  = 0.0

        # "Shuffle" our training data so each epoch sees it in a different
        # order.  Note that we generate a permutation of batch indices so
        # we don't actually rearrange the training data in memory and
        # dramatically slow things down.
        #
        # NOTE: This is shuffling batches and *not* the individual training
        #       samples.  This will not be sufficient if the training samples
        #       are highly correlated (e.g. they are generated in a sequence)!
        #
        permuted_batch_indices = np.random.permutation( number_batches )

        for batch_index in range( number_batches ):
            start_index = permuted_batch_indices[batch_index] * batch_size
            end_index   = start_index + batch_size

            # Get the next batch of droplets.
            inputs          = training_inputs[start_index:end_index, :]
            outputs         = training_outputs[start_index:end_index, :]
            times           = training_times[start_index:end_index]
            current_weights = weights[start_index:end_index]

            # Normalize the inputs and outputs to [-1, 1].
            normalized_inputs  = normalize_droplet_parameters( inputs )
            normalized_outputs = normalize_droplet_parameters( outputs )

            # XXX: Need to rethink how we handle time as an input.  This is annoying
            #      to have to stack a reshaped vector each time.
            normalized_inputs = np.hstack( (normalized_inputs,
                                            times.reshape( (batch_size, 1) )) )

            normalized_inputs  = torch.from_numpy( normalized_inputs ).to( device )
            normalized_outputs = torch.from_numpy( normalized_outputs ).to( device )

            # Reset the parameter gradients.
            optimizer.zero_grad()

            # Perform the forward pass.
            normalized_approximations = model( normalized_inputs )

            # Estimate the loss - using weights if needed.
            if criterion == weighted_mse_loss:
                loss = weighted_mse_loss( normalized_approximations, normalized_outputs, current_weights )
            else:
                loss = criterion( normalized_approximations, normalized_outputs )

            # Backwards pass and optimization.
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

            # Detect when we're at the end of the mini-batch and update our
            # loss.  This properly handles the case where we have exactly one
            # mini batch-worth of training data.
            if batch_index % MINI_BATCH_SIZE == (MINI_BATCH_SIZE - 1):
                running_loss /= MINI_BATCH_SIZE

                training_loss.append( running_loss )
                training_loss_history.append( running_loss )

                running_loss = 0.0
        else:
            # Handle the case where we don't have enough data to finish a
            # mini-batch.
            running_loss /= number_batches

            training_loss.append( running_loss )
            training_loss_history.append( running_loss )

            running_loss = 0.0

        # We finished all of the batches.  Adjust the learning rate before we
        # checkpoint so it can be loaded and training resumed without additional
        # preparation.
        for parameter_group in optimizer.param_groups:
            parameter_group["lr"] *= lr_scale

        # Checkpoint if requested.
        if checkpoint_prefix is not None:
            save_model_checkpoint( checkpoint_prefix,
                                   epoch_index,
                                   model,
                                   optimizer,
                                   criterion,
                                   training_loss_history )

        # Measure our validation loss if we have data.
        if validation_file is not None:
            with torch.no_grad():
                running_loss = 0.0

                # Compute the number of batches we have, possibly including a
                # partial last batch.
                number_validation_batches = (validation_inputs.shape[0] + VALIDATION_BATCH_SIZE - 1) // VALIDATION_BATCH_SIZE

                # Evaluate each of the batches and accumulate a running loss.
                for validation_batch_index in range( number_validation_batches ):
                    start_index = validation_batch_index * VALIDATION_BATCH_SIZE
                    if validation_batch_index == number_validation_batches - 1:
                        end_index = -1
                    else:
                        end_index = start_index + VALIDATION_BATCH_SIZE

                    validation_approximations = model( validation_inputs[start_index:end_index].to( device ) )

                    # Estimate the loss - using weights if needed.
                    if criterion == weighted_mse_loss:
                        loss = weighted_mse_loss( validation_approximations,
                                                  validation_outputs[start_index:end_index].to( device ),
                                                  current_weights )
                    else:
                        loss = criterion( validation_approximations,
                                          validation_outputs[start_index:end_index].to( device ) )

                    running_loss += loss.item()

            validation_loss = running_loss / number_validation_batches

            validation_loss_history.append( validation_loss )
        else:
            # Provide a placeholder for the validation loss but don't add it to
            # the history.
            validation_loss = np.nan

        # Execute each of the callbacks.
        for epoch_callback in epoch_callbacks:
            epoch_callback( model,
                            epoch_index,
                            optimizer,
                            training_loss,
                            validation_loss )

    return training_loss_history, validation_loss_history

def weighted_mse_loss( inputs, targets, weights ):
    """
    Calculates MSE error between inputs and targets with a weight for each difference.

    Takes 3 arguments:

      inputs      - Array of any size.
      targets     - Array with shape matching inputs.
      weights     - Array with same shape as inputs and targets with coefficients for
                    each difference.

    Returns 1 value:

      loss - The weighted MSE error with the same shape as inputs.

    """
    return (weights * ((inputs - targets)**2)).mean()
