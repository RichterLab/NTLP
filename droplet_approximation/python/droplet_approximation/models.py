import numpy as np
import torch
import torch.nn as nn

from .data import read_training_file
from .physics import dydt,\
                     get_parameter_ranges, \
                     normalize_droplet_parameters, \
                     scale_droplet_parameters, \
                     solve_ivp_float32_outputs

class ResidualNet( nn.Module ):
    """
    4-layer multi-layer perceptron (MLP) with ReLU activations.  This aims to
    balance parameter count vs computational efficiency so that inferencing with
    it is faster than Gauss-Newton iterative solvers.

    This residual network learns the delta between the input particle size and
    temperature and the outputs, given the provided background conditions.
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
        out = torch.relu( self.fc1( x ) )
        out = torch.relu( self.fc2( out ) )
        out = torch.relu( self.fc3( out ) )
        out = self.fc4( out )

        out += x[..., 0:2]

        return out

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

def do_inference( input_parameters, times, model, device ):
    """
    Estimates droplet parameters using a specified model.  Model evaluation is
    performed on the CPU

    Takes 4 arguments:

      input_parameters - NumPy array, sized number_droplets x 6, containing the
                         input parameters for number_droplets-many droplets.  These are
                         provided in their natural, physical ranges.
      times            - NumPy array, sized number_droplets, containing the integration
                         times to evaluate each droplet at.
      model            - PyTorch model to use.
      device           - Device string to perform the evaluation on.

    Returns 1 value:

      output_parameters - NumPy array, sized number_droplets x 2, containing the
                          estimated radius and temperature for each of the droplets
                          at the specified integraition times.  These are in their
                          natural physical ranges.

    """

    eval_model = model.to( device )
    eval_model.eval()

    normalized_inputs = np.hstack( (normalize_droplet_parameters( input_parameters ),
                                    times.reshape( (-1, 1) )) ).astype( "float32" )
    # XXX: is this to(device) necessary?
    normalized_outputs = eval_model( torch.from_numpy( normalized_inputs ).to( device ) ).to( device ).detach().numpy()

    return scale_droplet_parameters( normalized_outputs )

def do_iterative_bdf( input_parameters, times ):
    """
    Evaluates a particle trajectory along given background inputs with
    bdf. Requires all inputs to be sorted with respect to time.

    Evalutes iteratively, using the output time/radius for the (n-1)th time for
    the n-th time. Background parameters are tracked with `input_parameters[2:]`
    for each time step.

    Takes 2 arguments:

      input_parameters    - NumPy array, sized number_time_steps x 6, containing the
                            input parameters for a single particle in order by time.
                            These are provided in their natural, physical ranges.
      times               - NumPy array, shaped 1 x data length containing
                            the time at each step.

    Returns 1 value:

      output_parameters   - NumPy array, sized input_parameters.shape[0] x 2, containing the
                            estimated trajectory of a particle integrated
                            along its background parameters with BDF
                            in natural ranges. Does NOT calculate
                            output for the last column so that the lengths
                            of the input and output match

    """

    integration_times = np.diff( times )

    output_parameters       = np.zeros( (input_parameters.shape[0], 2), dtype=np.float32 )
    output_parameters[0, :] = input_parameters[0, :2]

    # Evaluate
    for time_index in range( 1, input_parameters.shape[0] ):
        output_parameters[time_index, :] = solve_ivp_float32_outputs( dydt,
                                                                      [0, integration_times[time_index - 1]],
                                                                      output_parameters[time_index - 1, :],
                                                                      method="BDF",
                                                                      t_eval=[integration_times[time_index-1]],
                                                                      args=(input_parameters[time_index-1, 2:],) ).y[:, 0]

    return output_parameters

def do_iterative_inference( input_parameters, times, model, device ):
    """
    Estimates a particle trajectory along given background inputs. Requires all
    inputs to be sorted with respect to time.

    Evaluates iteratively, using the output time/radius for the (n-1)th time for
    the n-th time. Background parameters are tracked with `input_parameters[2:]`
    for each time step.

    Takes 4 arguments:

      input_parameters    - Array, sized number_time_steps x 6, containing the
                            input parameters for a single particle in order by time.
                            These are provided in their natural, physical ranges.
      times               - Array, sized 1 x data length containing the time at
                            each step.
      model               - PyTorch model to use.
      device              - Device string to evaluate on.

    Returns 1 value:

      output_parameters   - NumPy array, sized input_parameters.shape[0] x 2,
                            containing the estimated trajectory of a particle
                            integrated along its background parameters with the
                            MLP in natural ranges.  Does NOT calculate output
                            for the last column so that the lengths of the input
                            and output match.

    """

    eval_model = model.to( device )
    eval_model.eval()

    normalized_data   = np.array( normalize_droplet_parameters( input_parameters ) )
    integration_times = np.diff( times )

    output_parameters       = np.zeros( (normalized_data.shape[0], 2), dtype=np.float32 )
    output_parameters[0, :] = normalized_data[0, :2]

    normalized_inputs = np.hstack( (output_parameters[0, :],
                                    normalized_data[0, 2:],
                                    [integration_times[0]] ) ).astype( "float32" )

    # Evaluate
    for time_index in range( 1,  normalized_data.shape[0] ):
        output_parameters[time_index, :] = eval_model( torch.from_numpy( normalized_inputs ).to( device ) ).detach().numpy()

        if time_index < normalized_data.shape[0] - 1:
            normalized_inputs = np.hstack( (output_parameters[time_index, :],
                                            normalized_data[time_index, 2:],
                                            [integration_times[time_index]] ) ).astype( "float32" )

    return scale_droplet_parameters( output_parameters )

def generate_fortran_module( model_name, model_state, output_path ):
    """
    Creates a Fortran 2003 module that allows use of the supplied model with a
    batch size of 1 during inference.

    NOTE: This currently expects the supplied model state to represent a 4-layer
          MLP with a specific set weights/biases names.

    Takes 3 arguments:

      model_name      - Name of the model to write in the generated module's comments
                        so as to identify where the weights came from.
      model_state     - PyTorch model state dictionary for the model to expose.
      output_path     - File to write to.  If this exists it is overwritten.

    Returns nothing.

    """

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

    def write_module_prolog( model_name, model_state, output_fp ):
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

        # Get the current particle data ranges.
        parameter_ranges = get_parameter_ranges()

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
    !         input_parameters = [radius, temperature, salt_mass, air_temperature, rh, rhoa, t_final]
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

    ! Indices into the input and output vectors.  The input vector uses all of
    ! the indices while the output vector only uses the first two.
    integer, parameter :: RADIUS_INDEX          = 1
    integer, parameter :: TEMPERATURE_INDEX     = 2
    integer, parameter :: SALT_MASS_INDEX       = 3
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
    real*4, parameter :: SALT_MASS_LOG_RANGE(2)   = [{salt_mass_start:.2f}, {salt_mass_end:.2f}]
    real*4, parameter :: AIR_TEMPERATURE_RANGE(2) = [{air_temperature_start:.1f}, {air_temperature_end:.1f}]
    real*4, parameter :: RH_RANGE(2)              = [{rh_start:.2f}, {rh_end:.2f}]
    real*4, parameter :: RHOA_RANGE(2)            = [{rhoa_start:.2f}, {rhoa_end:.2f}]

    real*4, parameter :: RADIUS_LOG_MEAN      = SUM( RADIUS_LOG_RANGE ) / 2
    real*4, parameter :: TEMPERATURE_MEAN     = SUM( TEMPERATURE_RANGE ) / 2
    real*4, parameter :: SALT_MASS_LOG_MEAN   = SUM( SALT_MASS_LOG_RANGE ) / 2
    real*4, parameter :: AIR_TEMPERATURE_MEAN = SUM( AIR_TEMPERATURE_RANGE ) / 2
    real*4, parameter :: RH_MEAN              = SUM( RH_RANGE ) / 2
    real*4, parameter :: RHOA_MEAN            = SUM( RHOA_RANGE ) / 2

    !
    ! NOTE: We have be careful about lines no longer than 132 characters so we
    !       eliminate whitespace in the expressions.
    !
    real*4, parameter :: RADIUS_LOG_WIDTH      = (RADIUS_LOG_RANGE(2)-RADIUS_LOG_RANGE(1))/2
    real*4, parameter :: TEMPERATURE_WIDTH     = (TEMPERATURE_RANGE(2)-TEMPERATURE_RANGE(1))/2
    real*4, parameter :: SALT_MASS_LOG_WIDTH   = (SALT_MASS_LOG_RANGE(2)-SALT_MASS_LOG_RANGE(1))/2
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
    creation_date=current_date_time_str,
    number_inputs=model_state["fc1.weight"].shape[1],
    number_outputs=model_state["fc4.weight"].shape[0],
    number_layer1_neurons=model_state["fc1.weight"].shape[0],
    number_layer2_neurons=model_state["fc2.weight"].shape[0],
    number_layer3_neurons=model_state["fc3.weight"].shape[0],
    radius_start=parameter_ranges["radius"][0],
    radius_end=parameter_ranges["radius"][1],
    temperature_start=parameter_ranges["temperature"][0],
    temperature_end=parameter_ranges["temperature"][1],
    salt_mass_start=parameter_ranges["salt_mass"][0],
    salt_mass_end=parameter_ranges["salt_mass"][1],
    air_temperature_start=parameter_ranges["air_temperature"][0],
    air_temperature_end=parameter_ranges["air_temperature"][1],
    rh_start=parameter_ranges["relative_humidity"][0],
    rh_end=parameter_ranges["relative_humidity"][1],
    rhoa_start=parameter_ranges["rhoa"][0],
    rhoa_end=parameter_ranges["rhoa"][1]
        )

        print( preamble_str, file=output_fp )

    def write_model_initializaiton_routine( model_state, output_fp, indentation_str ):
        """
        Writes the model's initialization routine to the supplied file handle with
        a specific amount of indentation.

        Takes 3 arguments:

          model_state     - PyTorch model state dictionary for the model to initialize.
          output_fp       - File handle to write to.
          indentation_str - String prepended to each line of the subroutine written.

        Returns nothing.

        """

        def write_model_biases( model_state, output_fp, indentation_str ):
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

        def write_model_weights( model_state, output_fp, indentation_str ):
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

        # Our initialization routine is inside of a "contains" block so we're
        # indented one level deeper than the caller.
        inner_indentation_str = "    {:s}".format( indentation_str )

        # Prolog and epilog of our initialization routine, along with extra
        # newlines so we survive the application of .splitlines() below.
        initialization_prolog_str = "subroutine initialize_model()\n\n    implicit none\n\n"
        initialization_epilog_str = "\nend subroutine initialize_model"

        # Write the subroutine's prolog.
        for initialization_line in initialization_prolog_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( initialization_line ) > 0 else "",
                initialization_line ),
                file=output_fp )

        # Write the weights and biases separated by an empty line.
        write_model_weights( model_state, output_fp, inner_indentation_str )
        print( "", file=output_fp )
        write_model_biases( model_state, output_fp, inner_indentation_str )

        # Write the subroutine's epilog.
        for initialization_line in initialization_epilog_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( initialization_line ) > 0 else "",
                initialization_line ),
                file=output_fp )

    def write_model_inference_routine( model_state, output_fp, indentation_str ):
        """
        Writes the model's inference routine to the supplied file handle with
        a specific amount of indentation.

        Takes 3 arguments:

          model_state     - Unused argument kept for a consistent signature.
          output_fp       - File handle to write to.
          indentation_str - String prepended to each line of the subroutine written.

        Returns nothing.

        """

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
    normalized_input(SALT_MASS_INDEX)       = (log10( input(SALT_MASS_INDEX) ) - SALT_MASS_LOG_MEAN) / SALT_MASS_LOG_WIDTH
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

    ! Scale the outputs to the expected ranges.
    output(RADIUS_INDEX)      = 10.0**(output(RADIUS_INDEX) * RADIUS_LOG_WIDTH + RADIUS_LOG_MEAN)
    output(TEMPERATURE_INDEX) = output(TEMPERATURE_INDEX) * TEMPERATURE_WIDTH + TEMPERATURE_MEAN

end subroutine estimate

"""

        for inference_line in inference_str.splitlines():
            print( "{:s}{:s}".format(
                indentation_str if len( inference_line ) > 0 else "",
                inference_line ),
                   file=output_fp )

    def write_module_epilog( model_state, output_fp ):
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

    with open( output_path, "w" ) as output_fp:
        # Write out the beginning of the module.  This includes the definitions for
        # the layers' weights and biases but not the method that initializes them.
        write_module_prolog( model_name, model_state, output_fp )

        # Write out both of the subroutines that comprise this model.
        write_model_initializaiton_routine( model_state, output_fp, indentation_str )
        write_model_inference_routine( model_state, output_fp, indentation_str )

        # Write out the end of the module.
        write_module_epilog( model_state, output_fp )

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

def train_model( model, criterion, optimizer, device, number_epochs, training_file ):
    """
    Trains the supplied model for one or more epochs using all of the droplet parameters
    in an on-disk training file.  The parameters are read into memory once and then
    randomly sampled each epoch.  Any weird parameters encountered in the training file
    are logged and

    Takes 6 arguments:

      model         - PyTorch model to optimize.
      criterion     - PyTorch loss object to use during optimization.
      optimizer     - PyTorch optimizer associated with model.
      device        - Device string indicating where the optimization is
                      being performed.
      number_epochs - Number of epochs to train model for.  All training
                      data in training_file will be seen by the model
                      this many times.
      training_file - Path to the file containing training data created by
                      create_training_file().

    Returns 1 value:

      loss_history - List of training losses, one per mini-batch during
                     the training process.

    """

    #
    # NOTE: This is inefficient and requires the entire training data set to reside in
    #       RAM.  We could read chunks of the file on demand but that would require
    #       a more sophisticated training loop that performs I/O in a separate thread
    #       while the optimization process executes.
    #
    #       TL;DR Find a big machine to train on.
    #
    input_parameters, output_parameters, integration_times = read_training_file( training_file )

    weights = np.reciprocal( integration_times )
    weights = np.stack( (weights, weights), axis=-1 )
    weights = torch.from_numpy( weights ).to( device )

    BATCH_SIZE      = 1024
    MINI_BATCH_SIZE = 1024

    #
    # NOTE: This ignores parameters if the last batch isn't complete.
    #
    NUMBER_BATCHES  = input_parameters.shape[0] // BATCH_SIZE

    # Track each mini-batch's training loss for analysis.
    loss_history = []

    batch_indices = np.arange( NUMBER_BATCHES )

    for epoch_index in range( number_epochs ):
        model.train()

        # Reset our training loss.
        running_loss = 0.0

        # "Shuffle" our training data so each epoch sees it in a different
        # order.  Note that we generate a permutation of batch indices so
        # we don't actually rearrange the training data in memory and
        # dramatically slow things down.
        permuted_batch_indices = np.random.permutation( batch_indices )

        for batch_index in range( NUMBER_BATCHES ):
            start_index = permuted_batch_indices[batch_index] * BATCH_SIZE
            end_index   = start_index + BATCH_SIZE

            # Get the next batch of droplets.
            inputs          = input_parameters[start_index:end_index, :]
            outputs         = output_parameters[start_index:end_index, :]
            times           = integration_times[start_index:end_index]
            current_weights = weights[start_index:end_index]

            # Normalize the inputs and outputs to [-1, 1].
            normalized_inputs  = normalize_droplet_parameters( inputs )
            normalized_outputs = normalize_droplet_parameters( outputs )

            # XXX: Need to rethink how we handle time as an input.  This is annoying
            #      to have to stack a reshaped vector each time.
            normalized_inputs = np.hstack( (normalized_inputs,
                                            times.reshape( (BATCH_SIZE, 1) )) )

            normalized_inputs  = torch.from_numpy( normalized_inputs ).to( device )
            normalized_outputs = torch.from_numpy( normalized_outputs ).to( device )

            # Reset the parameter gradients.
            optimizer.zero_grad()

            # Perform the forward pass.
            normalized_approximations = model( normalized_inputs )

            # Estimate the loss - using weights if needed
            if criterion == weighted_mse_loss:
                loss = weighted_mse_loss( normalized_approximations, normalized_outputs, current_weights )
            else:
                criterion( normalized_approximations, normalized_outputs )

            # Backwards pass and optimization.
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

            if batch_index > 0 and batch_index % MINI_BATCH_SIZE == 0:
                running_loss /= MINI_BATCH_SIZE

                print( "Epoch #{:d}, batch #{:d} loss: {:g}".format(
                    epoch_index + 1,
                    batch_index + 1,
                    running_loss ), flush=True )
                print( ((normalized_approximations - normalized_outputs)**2).mean( axis=0 ) )
                print( ((normalized_approximations - normalized_outputs)**2).max( axis=0 ) )
                loss_history.append( running_loss )

                # Break out of the batch loop if we ever encounter an input
                # that causes the loss to spike.  Training data generation
                # should produce good inputs and outputs though should something
                # slip through we want to immediately stop training so we can
                # understand what went wrong - there is no way to recover
                # from a loss spike that is O(10) when our target loss is O(1e-4).
                if running_loss > 10:
                    print( "Crazy loss!" )
                    print( inputs, outputs )
                    break

                running_loss = 0.0
        else:
            # We finished all of the batches.  Adjust the learning rate and
            # go to the next epoch.
            for parameter_group in optimizer.param_groups:
                parameter_group["lr"] /= 2

            continue

        # We didn't complete all of the batches in this epoch.
        break

    # Handle the case where we didn't have enough data to complete a mini-batch.
    if len( loss_history ) == 0:
        running_loss /= number_epochs * NUMBER_BATCHES
        loss_history.append( running_loss )

        print( "Epoch #{:d}, batch #{:d} loss: {:g}".format(
            number_epochs,
            NUMBER_BATCHES,
            running_loss ) )

    return loss_history

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
