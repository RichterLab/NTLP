import fcntl
import multiprocessing
import os
import warnings

import numpy as np
import pandas as pd

from .physics import TimeoutError, \
                     dydt, \
                     get_parameter_ranges, \
                     scale_droplet_parameters, \
                     timed_solve_ivp

def batch_convert_NTLP_traces_to_particle_files( trace_paths, particles_root, dirs_per_level, number_processes=0 ):
    """
    Converts a sequence of NTLP particle trace files into raw particle files and
    returns the path to a particle index file containing a sorted list of all
    unique particle identifiers seen during the conversion.  Conversion of
    multiple NTLP trace files is performed in parallel.

    Takes 4 arguments:

      trace_paths      - Sequence of NTLP particle trace paths to convert.
      particles_root   - Path to the top-level directory to write raw particle
                         files beneath.
      dirs_per_level   - Number of subdirectories per level in the particle files
                         directory hierarchy.
      number_processes - Optional number of processes to use for converting NTLP
                         particle trace files in parallel.  If omitted, defaults to
                         0 and one process per core will be utilized.

    Returns 2 values:

      particles_index_path - Path to the particles index file written.  This is beneath
                             particles_root.
      unique_particle_ids  - Sorted NumPy array of unique particle identifiers parsed
                             from all of paths in trace_paths.

    """

    particles_index_path = "{:s}/particles.index".format( particles_root )

    # Default to the number of cores on the system.
    #
    # NOTE: This may not be ideal if symmetric multi-threading (SMT, aka
    #       HyperThreading) is turned on.
    #
    if number_processes == 0:
        number_processes = multiprocessing.cpu_count()

    # Convert all of the particles in parallel.
    with multiprocessing.Pool( processes=number_processes ) as pool:
        args_list = [(trace_path, particles_root, dirs_per_level) for trace_path in trace_paths]

        particle_ids_list = pool.starmap( convert_NTLP_trace_to_particle_files,
                                          args_list )

    # Serialize each of the unique particle identifiers seen to the index.  This
    # provides discoverability to the raw particle files written above.
    unique_particle_ids = np.unique( np.hstack( particle_ids_list ) )
    unique_particle_ids.tofile( particles_index_path )

    return particles_index_path, unique_particle_ids

def batch_read_particles_data( particles_root, particle_ids, dirs_per_level, number_processes=0, quiet_flag=True ):
    """
    Reads one or more raw particle files and returns a particles DataFrame.
    Reading is performed in parallel.

    Takes 5 arguments:

      particles_root   - Path to the top-level directory to read raw particle
                         files from.
      particle_ids     - Sequence of particle identifiers to read.
      dirs_per_level   - Number of subdirectories per level in the particle
                         files directory hierarchy.
      number_processes - Optional number of processes to use for reading
                         particle trace files in parallel.  If omitted, defaults
                         to 0 and one process per core will be utilized.
      quiet_flag       - Optional flag specifying whether DataFrame creation
                         should be quiet or not.  If omitted, defaults to True
                         and creation does not generate outputs unless an error
                         occurs.

    Returns 1 value:

      particles_df - DataFrame with one row per partical identifier in particle_ids.
                     See read_particles_data() for details.

    """

    # Default to the number of cores on the system.
    #
    # NOTE: This may not be ideal if symmetric multi-threading (SMT, aka
    #       HyperThreading) is turned on.
    #
    if number_processes == 0:
        number_processes = multiprocessing.cpu_count()

    # Figure out how many particles to process per process.  The last process
    # gets fewer.
    number_particles             = len( particle_ids )
    number_particles_per_process = number_particles // number_processes

    # Create a list of ranges so we can easily construct the arguments to
    # the DataFrame parsing code below.
    particle_ranges = []
    for process_index in np.arange( number_processes ):
        start_index = number_particles_per_process * process_index
        end_index   = number_particles_per_process * (process_index + 1)

        # The last process gets whatever is remaining.
        if process_index == (number_processes - 1):
            end_index = -1

        particle_range = range( start_index, end_index )
        particle_ranges.append( particle_range )

    # Read each chunk of particles in parallel.
    with multiprocessing.Pool( processes=number_processes ) as pool:
        args_list = []
        for particle_range in particle_ranges:
            args_list.append( (particles_root,
                               particle_ids[particle_range],
                               dirs_per_level,
                               quiet_flag) )

        dfs = pool.starmap( read_particles_data, args_list )

    # Concatenate the resulting DataFrames.
    return pd.concat( dfs, axis=0 )

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

    # Get the current parameter ranges.
    parameter_ranges = get_parameter_ranges()

    #
    # NOTE: This isn't terribly efficient and assumes there is enough RAM to
    #       read, make a copy, and write back out.
    #
    input_parameters, output_parameters, times = read_training_file( file_name )

    bad_data_mask = ((output_parameters[:, 0] < 10.0**parameter_ranges["radius"][0] * (100 - 3) / 100) |
                     (output_parameters[:, 0] > 10.0**parameter_ranges["radius"][1] * (100 + 3) / 100) |
                     (output_parameters[:, 1] < parameter_ranges["temperature"][0] * (100 - 3) / 100) |
                     (output_parameters[:, 1] > parameter_ranges["temperature"][1] * (100 + 3) / 100))

    input_outputs = np.hstack( (input_parameters,
                                times.reshape( (-1, 1) ),
                                output_parameters) )

    # Overwrite the original file with only the good data.
    input_outputs[~bad_data_mask, :].tofile( file_name )

    # Return the counts of good and bad data so the caller knows what was done.
    return input_parameters.shape[0], bad_data_mask.sum()

def convert_NTLP_trace_to_particle_files( trace_path, particles_root, dirs_per_level ):
    """
    Converts an NTLP trace file into a one or more raw particle files.  NTLP
    trace files contain observations of any particles that were seen by the MPI
    rank that created the trace though it is easier to analyze particles when
    they are stored in files dedicated to individual particle observations.

    NOTE: This does not generate a unique particle identifier index!

    Takes 3 arguments:

      trace_path     - Path to the NTLP particle trace file to process.
      particles_root - Path to the top-level directory to write raw particle
                       files beneath.
      dirs_per_level - Number of subdirectories per level in the particle files
                       directory hierarchy.

    Returns 1 value:

      particle_ids - NumPy array of the unique particle identifiers in trace_path.

    """

    def make_particle_file_directories( particle_path ):
        """
        Creates all of the intermediate directories in the supplied raw particle
        file path.  This helper function

        Takes 1 argument:

          particle_path - Path to the raw particle file.

        Returns nothing.

        """

        # Create a directory for a (non-existent) file, but do not do anything
        # when provided a directory.
        if not os.path.isdir( particle_path ):
            parent_path = "/".join( particle_path.split( "/" )[:-1] )
            os.makedirs( parent_path, exist_ok=True )

        return

    # Each particle is comprised of 2x 32-bit integers and 7x 32-bit floats.
    RECORD_SIZE       = 9
    RECORD_SIZE_BYTES = RECORD_SIZE * 4

    # Indices of the record's fields.
    PARTICLE_ID_INDEX = 0

    # Read in a large block of records at once.
    NUMBER_RECORDS_PER_CHUNK = 1024 * 1024

    # List of NumPy 1D arrays containing particle identifier associated with
    # each of the trace file's chunks.
    particle_ids_list = []

    with open( trace_path, "rb" ) as trace_fp:

        # Iterate through the file reading a chunk of records at
        # a time.
        while True:
            raw_buffer = trace_fp.read( RECORD_SIZE_BYTES * NUMBER_RECORDS_PER_CHUNK )

            # We reached the end of the file in the previous iteration.
            if not raw_buffer:
                break

            # Report a partial record.  This gets dropped as we only process an
            # integral number of records.
            if len( raw_buffer ) % RECORD_SIZE_BYTES:
                print( "'{:s}' contains a partial, trailing record!".format( trace_path ) )

            number_records = len( raw_buffer ) // RECORD_SIZE_BYTES

            # Create a 32-bit float array of the records read along with a
            # 32-bit integer view so we can access both halves of the record.
            # We make a 32-bit integer view on this data below.
            chunk_array_fp32  = np.frombuffer( raw_buffer[:(number_records * RECORD_SIZE_BYTES)],
                                               dtype=np.float32 ).reshape( -1, RECORD_SIZE )
            chunk_array_int32 = chunk_array_fp32.view( np.int32 )

            # Loop through each of the unique particles and append their records
            # into a separate file.
            particle_ids = np.unique( chunk_array_int32[:, PARTICLE_ID_INDEX] )
            for particle_index, particle_id in enumerate( particle_ids ):
                particle_path = get_particle_file_path( particles_root,
                                                        particle_id,
                                                        dirs_per_level )
                make_particle_file_directories( particle_path )

                with open( particle_path, "ab" ) as particle_fp:
                    particle_mask = (chunk_array_int32[:, PARTICLE_ID_INDEX] == particle_id)

                    try:
                        # Grab the lock so we can write our particles.  This
                        # blocks until we acquire it.  We need this to ensure
                        # correct operation when processing multiple NTLP
                        # particle trace files in parallel as the same particle
                        # may occur in different MPI ranks' traces.
                        fcntl.flock( particle_fp.fileno(), fcntl.LOCK_EX )

                        # Append our records to the file.
                        chunk_array_fp32[particle_mask, :].tofile( particle_fp )

                    except IOError as e:
                        raise RuntimeError( "Failed to lock '{:s}' for write: {:s}".format(
                            particle_path,
                            str( e ) ) )

                    finally:
                        # Unlock the file for the next writer.
                        fcntl.flock( particle_fp.fileno(), fcntl.LOCK_UN )

            # Add these particle identifiers to the list of seen identifiers.
            particle_ids_list.append( particle_ids )

    # Return the unique identifiers seen in this trace.
    return np.unique( np.hstack( particle_ids_list ) )

def create_droplet_batch( number_droplets, number_evaluations=1, particle_temperature_distribution=None ):
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

      number_droplets                     - Number of droplets to generate parameters for.
      number_evaluations                  - Optional number of integration times to evaluate each
                                            parameter at so that a temporal window can be learned.  If
                                            omitted defaults to a single point in time per parameter.
      particle_temperature_distribution   - Optional tuple of two floats desribing the mean and location
                                            for a normal distribution. This distribution will sample the
                                            difference between particle temperature and air temperature.
                                            If none, generate air temperature as usual.


    Returns 5 values:

      random_inputs     - Array, sized number_droplets x 6, containing the droplets
                          radii, temperatures, salt mass, air temperature, relative
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

    # Get the current parameter ranges.
    parameter_ranges = get_parameter_ranges()

    random_inputs     = np.empty( (number_droplets * number_evaluations, 6), dtype=np.float32 )
    random_inputs[::number_evaluations, :] = scale_droplet_parameters( np.reshape( np.random.uniform( -1, 1, number_droplets*6 ),
                                                                       (number_droplets, 6) ).astype( "float32" ) )

    # Fix the particle temperature based on air temperature
    if particle_temperature_distribution is not None:
        random_inputs[::number_evaluations, 1] = (random_inputs[::number_evaluations, 3] +
                                                  np.random.normal( loc=particle_temperature_distribution,
                                                                    scale=particle_temperature_distribution,
                                                                    size=number_droplets ))

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
            y0         = random_inputs[droplet_index, :2]
            parameters = random_inputs[droplet_index, 2:]
            t_final    = integration_times[droplet_index]

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
                    if solution.y[0][0] < 10.0**parameter_ranges["radius"][0] * (100 - 3) / 100:
                        weird_outputs["radius_too_small"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False
                    elif solution.y[0][0] > 10.0**parameter_ranges["radius"][1] * (100 + 3) / 100:
                        weird_outputs["radius_too_large"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False

                    if solution.y[1][0] < parameter_ranges["temperature"][0] * (100 - 3) / 100:
                        weird_outputs["temperature_too_small"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                        good_parameters_flag = False
                    elif solution.y[1][0] > parameter_ranges["temperature"][1] * (100 + 3) / 100:
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
                    random_inputs[droplet_index, :][2] *= 1.0 + (0.0001 * np.random.choice( [-1, 1] ))
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

def create_training_file( file_name, number_droplets, weird_file_name=None, user_batch_size=None, particle_temperature_distribution=None ):
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
      3. Input salt mass, in kilograms
      4. Input air temperature, in Kelvin
      5. Input relative humidity, as a non-dimensional value with 100% humidity at 1.0
      6. Input rhoa, in non-dimensional units
      7. Input evaluation time, in seconds
      8. Output radius, in meters
      9. Output temperature, in Kelvin

    The file generated can be read via read_training_file().

    Takes 4 arguments:

      file_name                           - Path to the file to write the droplet parameters.  This
                                            is overwritten if it exists.
      number_droplets                     - The number of droplets to generate.
      weird_file_name                     - Optional path to write any weird parameters encountered
                                            during batch generation.  If specified an Excel spreadsheet
                                            file is written, otherwise weird parameters are silently
                                            ignored.
      user_batch_size                     - Optional batch size specifying the number of parameters
                                            to generate at once.  If omitted,
      particle_temperature_distribution   - Optional tuple of two floats desribing the mean and location
                                            for a normal distribution. This distribution will sample the
                                            difference between particle temperature and air temperature.
                                            If none, generate air temperature as usual. Passed to
                                            `create_droplet_batch`

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
             batch_weird_outputs) = create_droplet_batch( batch_size,
                                                          particle_temperature_distribution=particle_temperature_distribution )

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

def get_particle_file_path( particles_root, particle_id, dirs_per_level ):
    """
    Returns the path to the raw particles file associated with the specified
    particle identifier.  This simply returns the path and does not guarantee
    its existence.

    Raw particle files are stored in a multi-level hierarchy so as to distribute
    them across multiple directories instead of having thousands in a single
    directory.  The `dirs_per_level` argument specifies how many directories are
    used at each level.

    Takes 3 arguments:

      particles_root - Path to the top-level of the raw particle files directory.
      particle_id    - Integral particles identifier of the particle to return
                       its path.
      dirs_per_level - Integral number of directories per level in the hierarchy.

    Returns 1 value:

      particle_file_path - Path to the requested particles file.

    """

    # Convert the particle identifier into a pair of integers in the range of
    # [0, dirs_per_level).  We use these components as sub-directory names
    # beneath the supplied root so that individual particle files are
    # distributed across multiple directories instead of having all of them
    # reside in a single subdirectory.
    #
    # To ensure that we spread particles across all top-level subdirectories
    # (breadth-first) we use the second component as the first subdirectory
    # and the first component as the second subdirectory.  This works best
    # when working with sequential particle identifiers.
    level_1 = (particle_id // dirs_per_level) % dirs_per_level
    level_2 = particle_id % dirs_per_level

    return "{:s}/{:03d}/{:03d}/{:d}.raw".format(
        particles_root,
        level_1,
        level_2,
        particle_id )

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
                       "salt_mass",
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
                       "salt_mass",
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

def normalize_NTLP_data( df ):
    """
    Adds columns for normalized droplet parameters to a NTLP dataframe from
    `read_NTLP_data`.

    Takes 1 argument:

      df  - Pandas DataFrame from `read_NTLP_data`.

    Returns nothing
    """

    # Get the current parameter ranges.
    parameter_ranges = get_parameter_ranges()

    # Normalize inputs
    df["normalized input radius"]      = (np.log10(df["input radius"]) - np.mean( parameter_ranges["radius"] )) / (np.diff( parameter_ranges["radius"])/2)
    df["normalized input temperature"] = (df["input temperature"] - np.mean( parameter_ranges["temperature"] )) / (np.diff( parameter_ranges["temperature"] )/2)
    df["normalized salt mass"]         = (np.log10(df["salt mass"]) - np.mean( parameter_ranges["salt_mass"] )) / (np.diff( parameter_ranges["salt_mass"] )/2)
    df["normalized air temperature"]   = (df["air temperature"] - np.mean( parameter_ranges["air_temperature"] )) / (np.diff( parameter_ranges["air_temperature"] )/2)
    df["normalized relative humidity"] = (df["relative humidity"] - np.mean( parameter_ranges["relative_humidity"] )) / (np.diff( parameter_ranges["relative_humidity"] )/2)
    df["normalized air density"]       = (df["air density"] - np.mean( parameter_ranges["rhoa"] )) / (np.diff( parameter_ranges["rhoa"])/2)

    # Normalize outputs
    df["normalized output radius"]      = (np.log10(df["output radius"]) - np.mean( parameter_ranges["radius"] )) / (np.diff( parameter_ranges["radius"] )/2)
    df["normalized output temperature"] = (df["output temperature"] - np.mean( parameter_ranges["temperature"] )) / (np.diff( parameter_ranges["temperature"] )/2)

def read_NTLP_data( file_name ):
    """
    Reads all of the fixed-size binary records from the path specified and returns
    a pandas DataFrame containing the processed data.

    Takes 1 argument:

      file_name           - String, path to NTLP dump

    Returns 1 value:

      df  - Pandas DataFrame containing:
        particle id         - integer, equals 100 x particle id + processor rank
        be flag             - integer flag, 0 if backwards euler failed to produce
                              an output. 1 if it successfully produced an output.
        time                - float, simulation time at which this sample was taken
        input radius        - float, radius of the particle
        input temperatur    - float, temperature of the particle
        salt mass           - float, salt mass of the particle
        air temperature     - float, ambient air temperature around particle
        relative humidity   - float, relative humidty at particle's surface
        air density         - float, air density around the particle
        integration time    - float, time between given time step until next time step
                              if `be_flag==1` and the next sample has the same particle.
                              Otherwise, `0`.
        output radius       - float, the radius at the next time step if `be flag==1`
                              and the next sample has the same particle. Otherwise, `0`.
        output temperature  - float, the temperature at the next time step if `be flag==1`
                              and the next sample has the same particle. Otherwise, `0`.

    """

    record_data_type = np.dtype( [("particle id",       np.int32),
                                  ("be flag",           np.int32),
                                  ("time",              np.float32),
                                  ("input radius",      np.float32),
                                  ("input temperature", np.float32),
                                  ("salt mass",         np.float32),
                                  ("air temperature",   np.float32),
                                  ("relative humidity", np.float32),
                                  ("air density",       np.float32)] )

    inputs_outputs = np.fromfile( file_name, dtype=record_data_type )

    df = pd.DataFrame( data=inputs_outputs,
                       columns=record_data_type.names )

    df["processor"] = df["particle id"] % 100

    df.sort_values( by=["particle id", "time"], ascending=True, inplace=True )

    # Calculate outputs and dt, set to 0 if be failed or new particle.
    #
    # NOTE: Failed backward Euler is denoted as a non-zero value.
    #
    calculate_output_flags = (df["be flag"] == 0) * (df["particle id"].diff( periods=-1 ) == 0)

    df["integration time"]   = calculate_output_flags * -df["time"].diff( periods=-1 )
    df["output radius"]      = calculate_output_flags *  df["input radius"].shift( periods=-1 )
    df["output temperature"] = calculate_output_flags *  df["input temperature"].shift( periods=-1 )

    return df

def read_particles_data( particles_root, particle_ids, dirs_per_level, quiet_flag=True ):
    """
    Reads one or more raw particle files and creates a particle-centric
    DataFrame with one row per particle read.  Each row contains all of the
    observations stored as 1D NumPy arrays making it easy to operate on the
    entirety of individual particle's lifetimes.

    Takes 4 arguments:

      particles_root - Path to the top-level directory to read raw particle
                       files from.
      particle_ids   - Sequence of particle identifiers to read the associated
                       raw particle files.
      dirs_per_level - Number of subdirectories per level in the particle files
                       directory hierarchy.
      quiet_flag     - Optional flag specifying whether DataFrame creation
                       should be quiet or not.  If omitted, defaults to True and
                       creation does not generate outputs unless an error
                       occurs.

    Returns 1 value:

      particle_df - DataFrame with one row per particle identifier in particle_ids.
                    Contains the following columns:

                      "number observations":  Number of observations for the particle.
                                              This is the length of each of the 1D array
                                              columns.
                      "birth time":           Simulation time when this particle was
                                              created.
                      "death time":           Simulation time when this particle was
                                              destroyed.
                      "integration times":    1D NumPy array of the timestep size
                      "number be failures":
                      "be statuses":
                      "input radii":
                      "output radii":
                      "input temperatures":
                      "output temperatures":
                      "salt masses":
                      "air temperatures":
                      "relative humidities":
                      "air densities":

    """

    # Each particle is comprised of 2x 32-bit integers and 7x 32-bit floats.
    RECORD_SIZE = 9

    RECORD_BE_STATUS_INDEX         = 1
    RECORD_TIME_INDEX              = 2
    RECORD_RADIUS_INDEX            = 3
    RECORD_TEMPERATURE_INDEX       = 4
    RECORD_SALT_MASS_INDEX         = 5
    RECORD_AIR_TEMPERATURE_INDEX   = 6
    RECORD_RELATIVE_HUMIDITY_INDEX = 7
    RECORD_AIR_DENSITY_INDEX       = 8

    data_dict    = {"particle id": particle_ids}
    particles_df = pd.DataFrame( data_dict ).set_index( "particle id" )

    number_particles = len( particle_ids )

    if not quiet_flag:
        print( "Parsing {:d} particle{:s}".format(
            number_particles,
            "" if number_particles == 1 else "s"
        ) )

    zeros_int32   = pd.Series( np.zeros( number_particles, dtype=np.int32 ),
                               index=particle_ids )
    zeros_float32 = pd.Series( np.zeros( number_particles, dtype=np.float32 ),
                               index=particle_ids )

    #
    # NOTE: Order matters here!  Put the more useful information first so it's
    #       easily visible when inspecting the DataFrame.
    #
    # NOTE: We must initialize the scalar columns with a default value versus
    #       specifying a Series object with a specific type.  The latter doesn't
    #       work for integral data types as NaN isn't representable as a
    #       "missing value" type which means Series are changed to a np.float64
    #       data type to accommodate the initial NaNs.
    #
    particles_df["number observations"] = zeros_int32
    particles_df["birth time"]          = zeros_float32
    particles_df["death time"]          = zeros_float32
    particles_df["integration times"]   = pd.Series( dtype=object )
    particles_df["number be failures"]  = zeros_int32
    particles_df["be statuses"]         = pd.Series( dtype=object )
    particles_df["input radii"]         = pd.Series( dtype=object )
    particles_df["output radii"]        = pd.Series( dtype=object )
    particles_df["input temperatures"]  = pd.Series( dtype=object )
    particles_df["output temperatures"] = pd.Series( dtype=object )
    particles_df["salt masses"]         = pd.Series( dtype=object )
    particles_df["air temperatures"]    = pd.Series( dtype=object )
    particles_df["relative humidities"] = pd.Series( dtype=object )
    particles_df["air densities"]       = pd.Series( dtype=object )

    for particle_id in particle_ids:
        particle_path = get_particle_file_path( particles_root, particle_id, dirs_per_level )

        # Get each of the observations.
        observations_fp32  = np.fromfile( particle_path, dtype=np.float32 ).reshape( -1, RECORD_SIZE )

        # Sort the observations by time.
        #
        # NOTE: This is needed in case this particle moved across MPI ranks
        #       during the simulation and was handled by more than one rank.
        #       In that case the observations within each rank are sorted but we
        #       could get windows of observations out of order depending on how
        #       the resulting trace files were post-processed into raw particle files.
        #
        sorted_indices     = np.argsort( observations_fp32[:, RECORD_TIME_INDEX] )
        observations_fp32  = observations_fp32[sorted_indices, :]
        observations_int32 = observations_fp32.view( np.int32 )

        if not quiet_flag:
            print( "Processing {:d}".format( int( particle_id ) ) )

        # Total number of full observations for this particle.  The last observation
        # is ignored since we're combining the inputs of successive observations as
        # the previous observation's outputs, and we don't have it's result.
        particles_df.at[particle_id, "number observations"] = observations_fp32.shape[0] - 1

        # Backward Euler failures for this particle.
        particles_df.at[particle_id, "number be failures"]  = (observations_int32[:-1, RECORD_BE_STATUS_INDEX] > 0).sum()
        particles_df.at[particle_id, "be statuses"]         = observations_int32[:-1, RECORD_BE_STATUS_INDEX]

        # Simulation timeline.
        particles_df.at[particle_id, "birth time"]          = observations_fp32[0,  RECORD_TIME_INDEX]
        particles_df.at[particle_id, "death time"]          = observations_fp32[-1, RECORD_TIME_INDEX]
        particles_df.at[particle_id, "integration times"]   = np.diff( observations_fp32[:,  RECORD_TIME_INDEX] )

        # Observed radii.
        particles_df.at[particle_id, "input radii"]         = observations_fp32[:-1, RECORD_RADIUS_INDEX]
        particles_df.at[particle_id, "output radii"]        = observations_fp32[1:,  RECORD_RADIUS_INDEX]

        # Observed temperatures.
        particles_df.at[particle_id, "input temperatures"]  = observations_fp32[:-1, RECORD_TEMPERATURE_INDEX]
        particles_df.at[particle_id, "output temperatures"] = observations_fp32[1:,  RECORD_TEMPERATURE_INDEX]

        # Background parameters.
        particles_df.at[particle_id, "salt masses"]         = observations_fp32[:-1, RECORD_SALT_MASS_INDEX]
        particles_df.at[particle_id, "air temperatures"]    = observations_fp32[:-1, RECORD_AIR_TEMPERATURE_INDEX]
        particles_df.at[particle_id, "relative humidities"] = observations_fp32[:-1, RECORD_RELATIVE_HUMIDITY_INDEX]
        particles_df.at[particle_id, "air densities"]       = observations_fp32[:-1, RECORD_AIR_DENSITY_INDEX]

    return particles_df

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
