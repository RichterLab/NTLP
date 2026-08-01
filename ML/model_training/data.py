from enum import IntEnum
import fcntl
import multiprocessing
import os
import warnings

import numpy as np
import pandas as pd

from physics import BDF_TOLERANCE_ABSOLUTE, \
                    BDF_TOLERANCE_RELATIVE, \
                    TimeoutError, \
                    dydt, \
                    get_parameter_ranges, \
                    generate_random_droplets_sobol, \
                    scale_droplet_parameters, \
                    timed_solve_ivp

def create_droplet_batch( number_droplets, number_evaluations=1, solve_ivp_atol=None, solve_ivp_rtol=None, max_salinity=np.inf, sobol_index=0, sobol_seed=0 ):
    """
    Creates a batch of random droplets' input parameters, with t_final.  t_final is sampled
    from a slightly larger distribution than the anticipated use cases (spanning [1e-3, 1e1])
    so as to increase the performance at the edges of the range.

    Takes 7 arguments:

      number_droplets    - Number of droplets to generate parameters for.
      number_evaluations - Optional number of integration times to evaluate each
                           parameter at so that a temporal window can be learned.  If
                           omitted defaults to a single point in time per parameter.
      solve_ivp_atol     - Optional scalar or pair of floating point values specifying
                           the absolute tolerance used by solve_ivp().  If a scalar is
                           provided, then both the radius and temperature tolerances
                           are set to solve_ivp_atol.  Otherwise, the first value
                           is the radius tolerance and the second is the temperature
                           tolerance.  If omitted, defaults to physics.BDF_TOLERANCE_ABSOLUTE.
      solve_ivp_rtol     - Scalar floating point value specifying the relative
                           tolerance used by solve_ivp().  If omitted, defaults to
                           physics.BDF_TOLERANCE_RELATIVE.
      max_salinity       - Optional scalar floating point, in the range of [0,
                           np.inf) specifying the maximum salinity allowed for a
                           droplet.  Droplets whose salt solute mass are greater
                           than this salinity will have their mass adjusted to a
                           smaller value that satisfies this maximum.  Salinity
                           is defined as the ratio of salt solute divided by the
                           droplet's water volume times the dimensionless
                           density of fresh water.  If omitted, defaults to
                           np.inf which corresponds to wet salt and admits any
                           salt solute mass.
      sobol_index        - Optional integer determining how far to fast forward the sobol sequence.
      sobol_seed         - Optional integer to seed the sobol sequence for data generation.

    Returns 3 values:

      random_inputs     - Array, sized number_droplets x 6, containing the droplets
                          radii, temperatures, salt solute, air temperature, relative
                          humidity, and rhoa.
      random_outputs    - Array, sized number_droplets x 2, containing the droplets
                          radii and temperatures.
      integration_times - Array, sized number_droplets, containing the times corresponding
                          to the associated random_inputs and random_outputs.

    """

    if solve_ivp_atol is None:
        absolute_tolerance = BDF_TOLERANCE_ABSOLUTE
    else:
        absolute_tolerance = solve_ivp_atol

    if solve_ivp_rtol is None:
        relative_tolerance = BDF_TOLERANCE_RELATIVE
    else:
        relative_tolerance = solve_ivp_rtol

    # Promote run-time warnings to errors so we get an exception whenever the ODEs
    # are evaluated in problematic corners of the parameter space.  This not only
    # prevents them from showing up on standard error but also lets us log them and
    # ignore them so we only sample the "good" parts of the space.
    warnings.simplefilter( "error", RuntimeWarning )

    # Get the current parameter ranges.
    parameter_ranges = get_parameter_ranges()

    random_inputs                          = np.empty( (number_droplets * number_evaluations, 6),
                                                       dtype=np.float32 )
    random_inputs[::number_evaluations, :] = generate_random_droplets_sobol( number_droplets,
                                                                             sobol_index=sobol_index,
                                                                             sobol_seed=sobol_seed,
                                                                             max_salinity=max_salinity )

    # Duplicate each unique droplet parameter once for each evaluation.
    # This keeps them in parameter order.
    for droplet_index in np.arange( number_droplets ):
        start_index = droplet_index * number_evaluations
        end_index   = (droplet_index + 1) * number_evaluations

        random_inputs[start_index+1:end_index, :] = random_inputs[start_index, :]

    random_outputs = np.empty_like( random_inputs, shape=(number_droplets * number_evaluations, 2) )

    # Size of the time window to sample when we're generating multiple temporal
    # evaluations.  This is one order of magnitude.
    TIME_WINDOW_SIZE = 1

    # We generate data for a time window that is slightly larger than what we're interested
    # in (logspace( -3, 1 )) so we can learn the endpoints.
    TIME_RANGE = (10.0**parameter_ranges["time"][0], 10.0**parameter_ranges["time"][1])

    integration_times = np.empty_like( random_inputs, shape=(number_droplets * number_evaluations) )

    # Generate the starting point for each of the evaluation groups.  We
    # construct the other times in the group below.
    integration_times[::number_evaluations] = (10.0**np.random.uniform( np.log10( TIME_RANGE[0] ),
                                                                        np.log10( TIME_RANGE[1] ),
                                                                        number_droplets )).astype( "float32" )

    if number_evaluations > 1:
        for droplet_index in np.arange( number_droplets ):
            start_index = droplet_index * number_evaluations
            end_index   = (droplet_index + 1) * number_evaluations

            # Generate the remaining points in this group while taking care
            # to not go beyond the largest time we're allowed.
            starting_exponent = np.log10( integration_times[start_index] )
            integration_times[start_index:end_index] = 10.0**np.linspace( starting_exponent,
                                                                          min( starting_exponent + TIME_WINDOW_SIZE,
                                                                               np.log10( TIME_RANGE[1] ) ),
                                                                          number_evaluations )

    # XXX: Need to figure out if we can simply evaluate the ODEs at number_evaluations-many points
    #      and re-roll only the points that are outside the expected ranges.  Right now we evaluate
    #      the ODEs for every t_final even when we could evaluate it once for a list of multiple so
    #      that we discard all of the evaluations when one of them is problematic.
    for droplet_index in np.arange( number_droplets * number_evaluations ):

        # Number of times the inputs have been nudged because they resulted in a
        # timeout during the ODE solve.
        nudge_count = 0

        # Emulate a do/while loop so we always evaluate at least one parameter before
        # deciding whether to keep it or not.
        while True:
            y0         = random_inputs[droplet_index, :2]
            parameters = random_inputs[droplet_index, 2:]
            t_final    = integration_times[droplet_index]

            try:
                random_outputs[droplet_index, :] = timed_solve_ivp( dydt,
                                                                    [0, t_final],
                                                                    y0,
                                                                    method="BDF",
                                                                    t_eval=[t_final],
                                                                    args=(parameters,),
                                                                    atol=absolute_tolerance,
                                                                    rtol=relative_tolerance )

                # Did the solve fail (e.g. dydt's exponentials overflowed, BDF's
                # LU solve failed, etc)?
                if not np.any( np.isnan( random_outputs[droplet_index, :] ) ):
                    good_parameters_flag = True

                    # Check that we didn't get a physically impossible solution.
                    # We reroll the dice to replace them.
                    #
                    # NOTE: We don't strictly need to check for negative temperatures as
                    #       that will get covered in validate_output_parameters() but
                    #       we leave it here so it is easy to identify physically impossible
                    #       cases.
                    #
                    if random_outputs[droplet_index, 0] <= 0.0:
                        good_parameters_flag = False
                    elif random_outputs[droplet_index, 1] <= 0.0:
                        good_parameters_flag = False

                    # Check that we didn't get strange solutions that are outside of the
                    # expected ranges.  These will also be logged and replaced.
                    if random_outputs[droplet_index, 0] < 10.0**parameter_ranges["radius"][0] * (100 - 3) / 100:
                        good_parameters_flag = False
                    elif random_outputs[droplet_index, 0] > 10.0**parameter_ranges["radius"][1] * (100 + 3) / 100:
                        good_parameters_flag = False

                    if random_outputs[droplet_index, 1] < parameter_ranges["temperature"][0] * (100 - 3) / 100:
                        good_parameters_flag = False
                    elif random_outputs[droplet_index, 1] > parameter_ranges["temperature"][1] * (100 + 3) / 100:
                        good_parameters_flag = False

                    # Jump to the next droplet's parameters if we didn't detect a problem.
                    if good_parameters_flag:
                        break

            except TimeoutError as e:
                nudge_count += 1
                if nudge_count < 3:
                    # Adjust the salt solute by 0.01% and see if that gets past
                    # whatever numerical issue the ODE solver has encountered.
                    random_inputs[droplet_index, :][2] *= 1.0 + (0.0001 * np.random.choice( [-1, 1] ))
                    continue
                else:
                    # Fall through and try a completely different set of parameters
                    # if we couldn't quickly find a solution.
                    nudge_count = 0

            # We failed to create acceptable parameters.  Reroll the dice for
            # this droplet and try again.
            random_inputs[droplet_index, :] = generate_random_droplets( 1,
                                                                        max_salinity=max_salinity )

    return random_inputs, random_outputs, integration_times

def create_training_file( file_name, number_droplets, user_batch_size=None, quiet_flag=True, sobol_seed=0, file_index=0 ):
    """
    Generates random droplet parameters, both inputs and their corresponding
    ODE outputs, and writes them as fixed-size binary records to a file.  This
    makes for efficient access, both sequential and random, for training and
    analysis.

    The output file written is comprised of one or more 36-byte records, one
    per droplet.  Each record holds 9x 32-bit, floating point in host-byte
    order (typically little-endian):

      1. Input radius, in meters
      2. Input temperature, in Kelvin
      3. Input salt solute, in kilograms
      4. Input air temperature, in Kelvin
      5. Input relative humidity, as a non-dimensional value with 100% humidity at 1.0
      6. Input rhoa, in non-dimensional units
      7. Input evaluation time, in seconds
      8. Output radius, in meters
      9. Output temperature, in Kelvin

    The file generated can be read via read_training_file().

    Takes 6 arguments:

      file_name       - Path to the file to write the droplet parameters.  This
                        is overwritten if it exists.
      number_droplets - The number of droplets to generate.
      user_batch_size - Optional batch size specifying the number of parameters
                        to generate at once.  If omitted,
      quiet_flag      - Optional flag specifying whether execution should be quiet.
                        If omitted, defaults to True.
      sobol_seed      - Optional integer to seed the sobol sequence for data generation.
      file_index      - Optional integer to fast forward the sobol sequence based on the
                        index of the file being generated.
    """

    # Balance the time it takes to generate a single batch vs the efficiency of
    # writing it out.  Each droplet's parameters takes 36 bytes (9x 32-bit floats
    # comprised of the 7x input parameters and the 2x output parameters).
    if user_batch_size is not None:
        BATCH_SIZE = user_batch_size
    else:
        BATCH_SIZE = 1024 * 8

    # Number of batches to create including the last, partial batch when
    # the batch size does not evenly divide the number of droplets.
    number_batches = (number_droplets + BATCH_SIZE - 1) // BATCH_SIZE

    with open( file_name, "wb" ) as output_fp:
        for batch_index in range( number_batches ):

            if not quiet_flag:
                print( "Writing batch #{:d} (configurations {:d}-{:d})".format(
                    batch_index,
                    batch_index * BATCH_SIZE,
                    (batch_index + 1) * BATCH_SIZE - 1 ) )

            # Determine how many droplets to create in this batch.
            if batch_index != (number_batches - 1):
                batch_size = BATCH_SIZE
            else:
                batch_size = number_droplets - BATCH_SIZE*batch_index

            # Get the next batch of droplets.
            (inputs,
             outputs,
             times) = create_droplet_batch( batch_size, sobol_index=batch_index*BATCH_SIZE + (number_droplets*file_index), sobol_seed=sobol_seed )

            inputs_outputs = np.hstack( (inputs,
                                         times.reshape( (batch_size, 1) ),
                                         outputs) )

            # Serialize the array
            inputs_outputs.tofile( output_fp )

def read_training_file( file_name ):
    """
    Reads all of the fixed-size binary records from the path specified and returns
    NumPy arrays containing input parameters and output parameters.  Note that the
    input parameters contains the integration times.

    Takes 1 arguments:

      file_name - Path to the file to parse.

    Returns 2 values:

      inputs  - NumPy array, shaped number_droplets x 7, containing the input parameters.
      outputs - NumPy array, shaped number_droplets x 2, containing the output parameters.

    """

    inputs_outputs = np.fromfile( file_name, dtype=np.float32 ).reshape( (-1, 9) )
    inputs         = inputs_outputs[:, :7]
    outputs        = inputs_outputs[:, 7:]

    return inputs, outputs
