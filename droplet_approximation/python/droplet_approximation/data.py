import warnings

import numpy as np

from .physics import DROPLET_RADIUS_LOG_RANGE, \
                     DROPLET_TEMPERATURE_RANGE, \
                     DROPLET_SALINITY_LOG_RANGE, \
                     DROPLET_AIR_TEMPERATURE_RANGE, \
                     DROPLET_RELATIVE_HUMIDITY_RANGE, \
                     DROPLET_RHOA_RANGE, \
                     TimeoutError, \
                     dydt, \
                     scale_droplet_parameters, \
                     timed_solve_ivp

def create_droplet_batch( number_droplets, linear_time_flag=False, number_evaluations=1 ):
    """
    Creates a batch of random droplets' input parameters, with t_final.  t_final is sampled
    from a slightly larger distribution than the anticipated use cases (spanning [1e-3, 1e1])
    so as to increase the performance at the edges of the range.

    Weird parameters encountered, both as inputs to the ODEs and the outputs returned, are
    logged and replaced with valid inputs and outputs to ensure the batch is suitable
    for training.  Logging returns lists of parameters categorized by the reason they were
    filtered out.

    The following categories exist for weird input parameters:

      evaluation_warning - The parameters generated a floating point exception during
                           evaluation of the ODEs
      failed_solve       - The ODEs failed to converge as a small enough timestep could
                           not be identified
      invalid_inputs     - XXX: Greg can't remember what triggered this and scipy.solve_ivp()
                                can raise ValueError in a number of ways

    The following categories exist for weird output parameters:

      radius_nonpositive      - The ODEs generated a physically impossible, negative radius
      radius_too_small        - The ODEs generated a radius smaller than the expected range
      radius_too_large        - The ODEs generated a radius larger than the expected range
      temperature_nonpositive - The ODEs generated a physically impossible, negative temperature
      temperature_too_small   - The ODEs generated a temperature smaller than the expected range
      temperature_too_large   - The ODEs generated a temperature larger than the expected range

    NOTE: The use of linear_time_flag==True leads to a poorly sampled t_final parameter
          that results in poor model performance for DNS time scales (i.e. t_final in
          [1e-3, 1e-1]) unless trained on a *lot* of droplets.  The default leads to
          decent performance on both DNS time scales and LES time scales (i.e. t_final
          in [1e-1, 1e1]).

    Takes 3 arguments:

      number_droplets    - Number of droplets to generate parameters for.
      linear_time_flag   - Optional boolean specifying whether integration times should
                           be sampled uniformly through the time range or log spaced.
                           Defaults to False so that short term dynamics are captured.
      number_evaluations - Optional number of integration times to evaluate each
                           parameter at so that a temporal window can be learned.  If
                           omitted defaults to a single point in time per parameter.

    Returns 5 values:

      random_inputs     - Array, sized number_droplets x 6, containing the droplets
                          radii, temperatures, salinity, air temperature, relative
                          humidity, and rhoa.
      random_outputs    - Array, sized number_droplets x 2, containing the droplets
                          radii and temperatures.
      integration_times - Array, sized number_droplets, containing the times corresponding
                          to the associated random_inputs and random_outputs.
      weird_inputs      - Dictionary containing lists of weird input parameters.  Each
                          key represents a category of "weird".
      weird_outputs     - Dictionary containing lists of weird output parameters.  Each
                          key represents a category of "weird".

    """

    # Promote run-time warnings to errors so we get an exception whenever the ODEs
    # are evaluated in problematic corners of the parameter space.  This not only
    # prevents them from showing up on standard error but also lets us log them and
    # ignore them so we only sample the "good" parts of the space.
    warnings.simplefilter( "error", RuntimeWarning )

    random_inputs     = np.empty( (number_droplets * number_evaluations, 6), dtype=np.float32 )
    random_inputs[::number_evaluations, :] = scale_droplet_parameters( np.reshape( np.random.uniform( -1, 1, number_droplets*6 ),
                                                                       (number_droplets, 6) ).astype( "float32" ) )
    # Duplicate each unique droplet parameter once for each evaluation.
    # This keeps them in parameter order.
    for droplet_index in np.arange( number_droplets ):
        start_index = droplet_index * number_evaluations
        end_index   = (droplet_index + 1) * number_evaluations

        random_inputs[start_index+1:end_index, :] = random_inputs[start_index, :]

    random_outputs = np.empty_like( random_inputs, shape=(number_droplets * number_evaluations, 2) )

    # Size of the time window to sample when we're generating multiple temporal
    # evaluations.  In linear units when linear_time_flag==True (e.g. 1 second),
    # otherwise in log-space (e.g. 1 order of magnitude).
    TIME_WINDOW_SIZE = 1

    # We generate data for a time window that is slightly larger than what we're interested
    # in (logspace( -3, 1 )) so we can learn the endpoints.
    TIME_RANGE = (10.0**-3.2, 10.0**1.1)

    integration_times = np.empty_like( random_inputs, shape=(number_droplets * number_evaluations) )
    if linear_time_flag:
        # Generate the starting point for each of the evaluation groups.  We
        # construct the other times in the group below.
        integration_times[::number_evaluations] = np.random.uniform( TIME_RANGE[0], TIME_RANGE[1], number_droplets ).astype( "float32" )

        if number_evaluations > 1:
            for droplet_index in np.arange( number_droplets ):
                start_index = droplet_index * number_evaluations
                end_index   = (droplet_index + 1) * number_evaluations

                # Generate the remaining points in this group while taking care
                # to not go beyond the largest time we're allowed.
                integration_times[start_index:end_index] = np.linspace( integration_times[start_index],
                                                                        min( integration_times[start_index] + TIME_WINDOW_SIZE, TIME_RANGE[1] ),
                                                                        number_evaluations )

    else:
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
                                                                              min( starting_exponent + TIME_WINDOW_SIZE, np.log10( TIME_RANGE[1] ) ),
                                                                              number_evaluations )

    # We count the number of problematic parameters we encounter.  Additionally,
    # track the reason why we found the parameters problematic along with the
    # parameters themselves.
    weird_inputs    = { "evaluation_warning": [],
                        "failed_solve":       [],
                        "invalid_inputs":     [],
                        "timeout_inputs":     [] }
    weird_outputs   = { "radius_nonpositive":      [],
                        "radius_too_small":        [],
                        "radius_too_large":        [],
                        "temperature_nonpositive": [],
                        "temperature_too_small":   [],
                        "temperature_too_large":   [] }

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
            y0         = random_inputs[droplet_index, :2] # Radius and temperature.
            parameters = random_inputs[droplet_index, 2:] # Environmental variables.
            t_final    = integration_times[droplet_index] # Time to integrate.

            try:
                solution     = timed_solve_ivp( dydt, [0, t_final], y0, method="BDF", t_eval=[t_final], args=(parameters,) )

                if solution.success:
                    random_outputs[droplet_index, 0] = solution.y[0][0]
                    random_outputs[droplet_index, 1] = solution.y[1][0]

                    good_parameters_flag = True

                    # Check that we didn't get a physically impossible solution.  These
                    # will be logged and we'll reroll the dice to replace them.
                    #
                    # NOTE: We don't strictly need to check for negative temperatures as
                    #       that will get covered in validate_output_parameters() but
                    #       we leave it here so it is easy to identify physically impossible
                    #       cases.
                    #
                    if solution.y[0][0] <= 0.0:
                        weird_outputs["radius_nonpositive"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False
                    elif solution.y[1][0] <= 0.0:
                        weird_outputs["temperature_nonpositive"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False

                    # Check that we didn't get strange solutions that are outside of the
                    # expected ranges.  These will also be logged and replaced.
                    if solution.y[0][0] < 10.0**DROPLET_RADIUS_LOG_RANGE[0] * (100 - 3) / 100:
                        weird_outputs["radius_too_small"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False
                    elif solution.y[0][0] > 10.0**DROPLET_RADIUS_LOG_RANGE[1] * (100 + 3) / 100:
                        weird_outputs["radius_too_large"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False

                    if solution.y[1][0] < DROPLET_TEMPERATURE_RANGE[0] * (100 - 3) / 100:
                        weird_outputs["temperature_too_small"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False
                    elif solution.y[1][0] > DROPLET_TEMPERATURE_RANGE[1] * (100 + 3) / 100:
                        weird_outputs["temperature_too_large"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False

                    # Jump to the next droplet's parameters if we didn't detect a problem.
                    if good_parameters_flag:
                        break
                else:
                    # Record the cases that fail to converge.
                    weird_inputs["failed_solve"].append( [y0, parameters, t_final, solution.message] )

            except RuntimeWarning as e:
                weird_inputs["evaluation_warning"].append( [y0, parameters, t_final, str( e )] )
            except TimeoutError as e:
                weird_inputs["timeout_inputs"].append( [y0, parameters, t_final, str( e )] )
                nudge_count += 1
                if nudge_count < 3:
                    # Adjust the salt content by 0.01% and see if that gets past
                    # whatever numerical issue the ODE solver has encountered.
                    random_inputs[droplet_index, :][2] *= 1.0 + (0.0001*np.random.choice( [-1, 1] ))
                    continue
                else:
                    # Fall through and try a completely different set of parameters
                    # if we couldn't quickly find a solution.
                    nudge_count = 0
            except ValueError as e:
                weird_inputs["invalid_inputs"].append( [y0, parameters, t_final, str( e )] )

            # We failed to create acceptable parameters.  Reroll the dice for
            # this droplet and try again.
            random_inputs[droplet_index, :] = scale_droplet_parameters(
                np.random.uniform( -1, 1, 6 ).astype( "float32" ) )

    # XXX: Warn users if there were strange droplets?  No simple way to see if any
    #      of the lists associated with weird_*'s keys are non-empty.

    return random_inputs, random_outputs, integration_times, weird_inputs, weird_outputs

def merge_weird_parameters( parameters_1, parameters_2 ):
    """
    Merge two dictionaries of lists into a separate copy containing the
    concatenated lists of both.  All list entries of the first dictionary
    come before the first entry of the second dictionary when the first
    has a non-empty list.

    Takes 2 arguments:

      parameters_1 - 1st dictionary of lists to merge.
      parameters_2 - 2nd dictionary of lists to merge.

    Returns 1 value:

      merged_parameers - Dictionary containing the merged lists.

    """

    merged_parameters = parameters_1.copy()

    for type_name, weird_things in parameters_2.items():
        if type_name in merged_parameters:
            merged_parameters[type_name].extend( weird_things )
        else:
            merged_parameters[type_name] = weird_things

    return merged_parameters

def write_weird_parameters_to_spreadsheet( file_name, weird_inputs, weird_outputs ):
    """
    Creates an Excel file containing one spreadsheet per type of weird inputs or
    outputs.

    See create_droplet_batch() for details on weird droplet inputs and outputs.

    Takes 3 arguments:

      file_name     - Path to the spreadsheet to create.  If it exists it is
                      overwritten.
      weird_inputs  - Dictionary of lists representing weird input droplet
                      parameters.
      weird_outputs - Dictionary of lists representing weird output droplet
                      parameters.

    Returns nothing.

    """

    # XXX: move this
    import pandas as pd

    # Open the file for writing, overwriting an existing file if necessary.
    with pd.ExcelWriter( file_name ) as writer:

        # Create a DataFrame for each category of input parameters weirdness.
        for type_name, weird_things in weird_inputs.items():
            sheet_name = "inputs-{:s}".format( type_name )

            # Columns in the DataFrame.
            columns = ["initial_radius",
                       "initial_temperature",
                       "salinity",
                       "air_temperature",
                       "relative_humidity",
                       "rhoa",
                       "t_final",
                       "error"]

            # Convert our list of lists of NumPy arrays to lists of lists
            # so spreadsheet cells contain a single value.
            weird_things_listified = []
            for weird_thing in weird_things:
                weird_thing_listified = []

                weird_thing_listified.extend( weird_thing[0].tolist() )
                weird_thing_listified.extend( weird_thing[1].tolist() )
                weird_thing_listified.extend( weird_thing[2:] )

                weird_things_listified.append( weird_thing_listified )

            # Build a DataFrame and write it as a new sheet.
            filtered_df = pd.DataFrame( weird_things_listified, columns=columns )
            filtered_df.to_excel( writer, sheet_name=sheet_name, index=False )

        # Create a DataFrame for each category of output parameters weirdness.
        for type_name, weird_things in weird_outputs.items():
            sheet_name = "outputs-{:s}".format( type_name )

            # Columns in the DataFrame.
            columns = ["initial_radius",
                       "initial_temperature",
                       "salinity",
                       "air_temperature",
                       "relative_humidity",
                       "rhoa",
                       "t_final",
                       "radius",
                       "temperature"]

            # Convert our list of lists of NumPy arrays to lists of lists
            # so spreadsheet cells contain a single value.
            weird_things_listified = []
            for weird_thing in weird_things:
                weird_thing_listified = []

                weird_thing_listified.extend( weird_thing[0].tolist() )
                weird_thing_listified.extend( weird_thing[1].tolist() )
                weird_thing_listified.append( weird_thing[2] )
                weird_thing_listified.extend( weird_thing[3].tolist() )

                weird_things_listified.append( weird_thing_listified )

            # Build a DataFrame and write it as a new sheet.
            filtered_df = pd.DataFrame( weird_things_listified, columns=columns )
            filtered_df.to_excel( writer, sheet_name=sheet_name, index=False )

def create_training_file( file_name, number_droplets, weird_file_name=None, user_batch_size=None ):
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
      3. Input salinity, aka the mass of disolved salt, in kilograms
      4. Input air temperature, in Kelvin
      5. Input relative humidity, as a non-dimensional value with 100% humidity at 1.0
      6. Input rhoa, in non-dimensional units
      7. Input evaluation time, in seconds
      8. Output radius, in meters
      9. Output temperature, in Kelvin

    The file generated can be read via read_training_file().

    Takes 4 arguments:

      file_name       - Path to the file to write the droplet parameters.  This
                        is overwritten if it exists.
      number_droplets - The number of droplets to generate.
      weird_file_name - Optional path to write any weird parameters encountered
                        during batch generation.  If specified an Excel spreadsheet
                        file is written, otherwise weird parameters are silently
                        ignored.
      user_batch_size - Optional batch size specifying the number of parameters
                        to generate at once.  If omitted, defaults to a "small"
                        number of parameters that balances memory footprint,
                        time to generate, and file write size.  This does not
                        need to evenly divide number_droplets.

    Returns nothing.

    """

    # Balance the time it takes to generate a single batch vs the efficiency of
    # writing it out.  Each droplet's parameters takes 36 bytes (9x 32-bit floats
    # comprised of the 7x input parameters and the 2x output parameters).
    if user_batch_size is not None:
        BATCH_SIZE = user_batch_size
    else:
        BATCH_SIZE = 1024 * 10

    # Number of batches to create including the last, partial batch when
    # the batch size does not evenly divide the number of droplets.
    number_batches = (number_droplets + BATCH_SIZE - 1) // BATCH_SIZE

    # Dictionaries for tracking weird parameters.
    weird_inputs  = {}
    weird_outputs = {}

    with open( file_name, "wb" ) as output_fp:
        for batch_index in range( number_batches ):

            print( "Writing batch #{:d} (configurations {:d}-{:d})".format(
                batch_index,
                batch_index * BATCH_SIZE,
                (batch_index + 1) * BATCH_SIZE - 1 ) )

            # Determine how many droplets to create in this batch.
            if batch_index != (number_batches - 1):
                batch_size = BATCH_SIZE
            else:
                batch_size = number_droplets % BATCH_SIZE

            # Get the next batch of droplets.
            (inputs,
             outputs,
             times,
             batch_weird_inputs,
             batch_weird_outputs) = create_droplet_batch( batch_size )

            # Track the weirdness for post-mortem analysis.
            weird_inputs  = merge_weird_parameters( weird_inputs, batch_weird_inputs )
            weird_outputs = merge_weird_parameters( weird_outputs, batch_weird_outputs )

            inputs_outputs = np.hstack( (inputs,
                                         times.reshape( (batch_size, 1) ),
                                         outputs) )

            # Serialize the array
            inputs_outputs.tofile( output_fp )

    if weird_file_name is not None:
        write_weird_parameters_to_spreadsheet( weird_file_name,
                                               weird_inputs,
                                               weird_outputs )

def read_training_file( file_name ):
    """
    Reads all of the fixed-size binary records from the path specified and returns
    NumPy arrays containing input parameters, output parameters, and integration
    times.

    Takes 1 arguments:

      file_name - Path to the file to parse.

    Returns 3 values:

      inputs  - NumPy array, shaped number_droplets x 6, containing the input parameters.
      outputs - NumPy array, shaped number_droplets x 2, containing the output parameters.
      times   - NumPy array, shaped number_droplets x 1, containing the integration times.

    """

    inputs_outputs = np.fromfile( file_name, dtype=np.float32 ).reshape( (-1, 9) )
    inputs         = inputs_outputs[:, :6]
    times          = inputs_outputs[:, 6]
    outputs        = inputs_outputs[:, 7:]

    return inputs, outputs, times

def clean_training_data( file_name ):
    """
    Removes invalid droplet parameters found in the supplied file and overwrites it
    so it only contains valid parameters.  Currently this only filters output
    parameters that are outside of the expected ranges by more than 3% on either end.

    NOTE: This function isn't normally needed, but there may be times where invalid
          outputs have been written but the vast majority are good and need to be
          salvaged.

    Takes 1 argument:

      file_name - Path to the file containing droplet parameters written by
                  create_training_file().  This file is overwritten with the
                  filtered droplet parameters.

    Returns 2 values:

      number_parameters     - Number of droplet parameters in the original file.
      number_bad_parameters - Number of droplet parameters that were filtered out.

    """

    #
    # NOTE: This isn't terribly efficient and assumes there is enough RAM to
    #       read, make a copy, and write back out.
    #
    input_parameters, output_parameters, times = read_training_file( file_name )

    bad_data_mask = ((output_parameters[:, 0] < 10.0**DROPLET_RADIUS_LOG_RANGE[0] * (100 - 3) / 100) |
                     (output_parameters[:, 0] > 10.0**DROPLET_RADIUS_LOG_RANGE[1] * (100 + 3) / 100) |
                     (output_parameters[:, 1] < DROPLET_TEMPERATURE_RANGE[0] * (100 - 3) / 100) |
                     (output_parameters[:, 1] > DROPLET_TEMPERATURE_RANGE[1] * (100 + 3) / 100))

    input_outputs = np.hstack( (input_parameters,
                                times.reshape( (-1, 1) ),
                                output_parameters) )

    # Overwrite the original file with only the good data.
    input_outputs[~bad_data_mask, :].tofile( file_name )

    # Return the counts of good and bad data so the caller knows what was done.
    return input_parameters.shape[0], bad_data_mask.sum()
