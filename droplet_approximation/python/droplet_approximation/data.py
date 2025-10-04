from enum import IntEnum
import fcntl
import multiprocessing
import os
import warnings

import numpy as np
import pandas as pd

from .physics import BDF_TOLERANCE_ABSOLUTE, \
                     BDF_TOLERANCE_RELATIVE, \
                     TimeoutError, \
                     dydt, \
                     get_parameter_ranges, \
                     scale_droplet_parameters, \
                     timed_solve_ivp

# "Evaluation" tag for backward Euler (BE) so it may be handled programmatically
# along with other evaluations.
BE_TAG_NAME = "be"

# Maximum integration gap size, in simulation seconds, to consider.  This is
# chosen such that it is larger than the largest LES-based timestep that
# could occur in NTLP.  Differences in observation times larger than one
# maximum timestep are are due to missing observations.
MAXIMUM_DT_GAP_SIZE = 10.0

class BEStatus( IntEnum ):
    """
    Container for constants describing the backward Euler (BE) status flags
    recorded in NTLP particle traces.
    """

    # BE was skipped due to the particle having a negative qinf value.  BE has
    # failed when this occurs (since BE could not complete successfully).
    NEG_QINF           = 1

    # The particle's radius was in equilibrium prior to BE.
    EQUILIBRIUM_RADIUS = 2

    # Gauss-Newton failed to converge, forcing at least one iteration of
    # Levenberg-Marquardt to execute.
    GN_FAILED          = 4

    # Levenberg-Marquardt failed to converge.  BE has failed when this occurred.
    LM_FAILED          = 8

    # Levenberg-Marquardt incorrectly converged and produced an invalid value.
    # This includes NaN and negative values for either radius or temperature,
    # as well as large radii.
    INVALID_OUTPUT     = 16

    # Flags locating the particle's position within the vertical domain
    # according to its input height.  One of the four flags is set based on
    # which quarter of the domain the particle was located (Q1 in [0, 25),
    # Q2 in [25, 50), Q3 in [50, 75), and Q4 in [75, 100]).
    HEIGHT_Q1          = 32
    HEIGHT_Q2          = 64
    HEIGHT_Q3          = 128
    HEIGHT_Q4          = 256

    # BE has failed when it could not successfully update the particle's radius
    # and temperature to a valid value.  Note that Levenberg-Marquardt runs if
    # Gauss-Newton fails to converge, so we only consider its failure as
    # critical.
    FAILED_UPDATE      = (NEG_QINF | LM_FAILED | INVALID_OUTPUT)

    def get_failure_mask( array, bit_mask=FAILED_UPDATE ):
        """
        Gets the Boolean array representing a bit mask applied to each element
        of an array.

        Takes 2 arguments:

          array    - Array of BE statues to locate failures in.
          bit_mask - Optional integral value, comprised of the logical OR of the
                     class' constants, to construct boolean_mask with.  If
                     omitted, defaults to FAILED_UPDATE.

        Returns 1 value:

          boolean_mask - Array of Boolean values indicating where the logical
                         AND of array and bit_mask was non-zero.

        """

        return (array & bit_mask) > 0

class ParticleRecord( IntEnum ):
    """
    Container for constants describing the format of a raw particle record.
    Provides the size of the record as well as named indices into the record
    itself.
    """

    # Indices in a raw particle record.  This is what is traced from NTLP; raw
    # particle files are comprised of records with the same particle identifier.
    PARTICLE_ID             = 0
    BE_STATUS_INDEX         = 1
    TIME_INDEX              = 2
    RADIUS_INDEX            = 3
    TEMPERATURE_INDEX       = 4
    SALT_SOLUTE_INDEX       = 5
    AIR_TEMPERATURE_INDEX   = 6
    RELATIVE_HUMIDITY_INDEX = 7
    AIR_DENSITY_INDEX       = 8

    # Each particle is comprised of 2x 32-bit integers and 7x 32-bit floats.
    SIZE       = AIR_DENSITY_INDEX + 1
    SIZE_BYTES = SIZE * 4

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

    particles_index_path = get_particles_index_path( particles_root )

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
    unique_particle_ids = write_particles_index( particles_index_path,
                                                 particle_ids_list )

    return particles_index_path, unique_particle_ids

def batch_read_particles_data( particles_root, particle_ids, dirs_per_level, number_processes=0, **kwargs ):
    """
    Reads one or more raw particle files and returns a particles DataFrame.
    Reading is performed in parallel.

    Takes 4 arguments:

      particles_root   - Path to the top-level directory to read raw particle
                         files from.
      particle_ids     - Sequence of particle identifiers to read.
      dirs_per_level   - Number of subdirectories per level in the particle
                         files directory hierarchy.
      number_processes - Optional number of processes to use for reading
                         particle trace files in parallel.  If omitted, defaults
                         to 0 and one process per core will be utilized.
      kwargs           - Optional dictionary of arguments to pass to
                         read_particles_data().

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
            end_index = number_particles

        particle_range = range( start_index, end_index )
        particle_ranges.append( particle_range )

    # Read each chunk of particles in parallel.
    with multiprocessing.Pool( processes=number_processes ) as pool:
        args_list = []
        for particle_range in particle_ranges:
            #
            # NOTE: We're passing the keyword arguments dictionary as a
            #       positional parameter to the internal wrapper!  This works
            #       around the fact that 1) .starmap() doesn't handle keyword
            #       arguments and 2) we can't run a lambda function (because it
            #       can't be pickled).
            #
            args_list.append( (particles_root,
                               particle_ids[particle_range],
                               dirs_per_level,
                               kwargs) )

        dfs = pool.starmap( _read_particles_data_wrapper, args_list )

    # Concatenate the resulting DataFrames.
    return pd.concat( dfs, axis=0 )

def be_status_to_z_domain_quartile( be_statuses ):
    """
    Converts backward Euler (BE) status flags into a the Z domain's quantile number
    where the evaluation started in.

    Takes 1 argument:

      be_statuses - NumPy array containing the BE status flags to convert.

    Returns 1 value:

      z_domain_quartile - NumPy array, shaped the same as be_statuses, containing the
                          Z domain's encoded quartile number.  Will be in the range of [1, 4]
                          representing the [0%, 25%), [25%, 50%), [50%, 75%), [75%, 100]%
                          quartiles.

    """

    # Construct a bitmask to pull out the Z domain flags.
    z_domain_bitmask = (BEStatus.HEIGHT_Q1 |
                        BEStatus.HEIGHT_Q2 |
                        BEStatus.HEIGHT_Q3 |
                        BEStatus.HEIGHT_Q4)

    quantiles = np.array( (be_statuses & z_domain_bitmask), dtype=np.int32 )

    # Remap the individual flags into the range of [1, 4] specifying the
    # quantile number the evaluation started in.
    for quantile_index in range( quantiles.shape[0] ):
        if quantiles[quantile_index] == BEStatus.HEIGHT_Q1:
            quantiles[quantile_index] = 1
        elif quantiles[quantile_index] == BEStatus.HEIGHT_Q2:
            quantiles[quantile_index] = 2
        elif quantiles[quantile_index] == BEStatus.HEIGHT_Q3:
            quantiles[quantile_index] = 3
        else:
            quantiles[quantile_index] = 4

    return quantiles

def be_success_mask( radius_data ):
    """
    Identifies trials where BE succeeded (True) and be failed (False).  This
    mask checks for trials that precede a lack of change in radius.  No change
    in radius corresponds to a time step with a failure in BE.

    Takes 1 argument:

      radius_data     - 1D NumPy array, length `number of time steps`,
                        containing the radius data for a particle sorted
                        by time.

    Returns 1 value:

      be_success_mask - 1D NumPy boolean array, of same length as radius_data,
                        containing True where BE suceeded, False where it
                        failed.

    """

    mask = radius_data[:-1] != radius_data[1:]

    #
    # NOTE: This assumes that BE always succeeds at the end!
    #
    return np.append( mask, True )

def convert_NTLP_trace_array_to_particle_files( trace_array, particles_root, dirs_per_level ):
    """
    Converts an array of NTLP trace observations into a one or more raw particle
    files.  The trace array contain observations of any particles that were seen
    by the MPI rank that created the trace though it is easier to analyze
    particles when they are stored in files dedicated to individual particle
    observations.

    NOTE: This does not generate a unique particle identifier index!

    Takes 3 arguments:

      trace_array    - NumPy array, shaped (number_observations,
                       ParticleRecord.SIZE), containing particle observations as
                       32-bit floating point values.
      particles_root - Path to the top-level directory to write raw particle
                       files beneath.
      dirs_per_level - Number of subdirectories per level in the particle files
                       directory hierarchy.

    Returns 1 value:

      particle_ids - NumPy array of the unique particle identifiers in trace_path.

    """

    def _make_particle_file_directories( particle_path ):
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

    # Create a 32-bit integer array of the records provided.  We need this to
    # determine which particles are present.
    trace_array_int32 = trace_array.view( np.int32 )

    # Loop through each of the unique particles and append their records
    # into a separate file.
    particle_ids = np.unique( trace_array_int32[:, ParticleRecord.PARTICLE_ID] )
    for particle_index, particle_id in enumerate( particle_ids ):
        particle_path = get_particle_file_path( particles_root,
                                                particle_id,
                                                dirs_per_level )
        _make_particle_file_directories( particle_path )

        with open( particle_path, "ab" ) as particle_fp:
            particle_mask = (trace_array_int32[:, ParticleRecord.PARTICLE_ID] == particle_id)

            try:
                # Grab the lock so we can write our particles.  This
                # blocks until we acquire it.  We need this to ensure
                # correct operation when processing multiple NTLP
                # particle trace files in parallel as the same particle
                # may occur in different MPI ranks' traces.
                fcntl.flock( particle_fp.fileno(), fcntl.LOCK_EX )

                # Append our records to the file.
                trace_array[particle_mask, :].tofile( particle_fp )

            except IOError as e:
                raise RuntimeError( "Failed to lock '{:s}' for write: {:s}".format(
                    particle_path,
                    str( e ) ) )

            finally:
                # Unlock the file for the next writer.
                fcntl.flock( particle_fp.fileno(), fcntl.LOCK_UN )

    return particle_ids

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

    # Read in a large block of records at once.
    NUMBER_RECORDS_PER_CHUNK = 1024 * 1024

    # List of NumPy 1D arrays containing particle identifier associated with
    # each of the trace file's chunks.
    particle_ids_list = []

    with open( trace_path, "rb" ) as trace_fp:

        # Iterate through the file reading a chunk of records at
        # a time.
        while True:
            raw_buffer = trace_fp.read( ParticleRecord.SIZE_BYTES * NUMBER_RECORDS_PER_CHUNK )

            # We reached the end of the file in the previous iteration.
            if not raw_buffer:
                break

            # Report a partial record.  This gets dropped as we only process an
            # integral number of records.
            if len( raw_buffer ) % ParticleRecord.SIZE_BYTES:
                print( "'{:s}' contains a partial, trailing record!".format( trace_path ) )

            number_records = len( raw_buffer ) // ParticleRecord.SIZE_BYTES

            # Create a NumPy array from our bytes, extract the particle
            # observations, and write them to disk.
            chunk_array_fp32 = np.frombuffer( raw_buffer[:(number_records * ParticleRecord.SIZE_BYTES)],
                                              dtype=np.float32 ).reshape( -1, ParticleRecord.SIZE )
            particle_ids     = convert_NTLP_trace_array_to_particle_files( chunk_array_fp32,
                                                                           particles_root,
                                                                           dirs_per_level )

            # Add these particle identifiers to the list of seen identifiers.
            particle_ids_list.append( particle_ids )

    # Return the unique identifiers seen in this trace.
    return np.unique( np.hstack( particle_ids_list ) )

def create_droplet_batch( number_droplets, number_evaluations=1, solve_ivp_atol=None, solve_ivp_rtol=None ):
    """
    Creates a batch of random droplets' input parameters, with t_final.  t_final is sampled
    from a slightly larger distribution than the anticipated use cases (spanning [1e-3, 1e1])
    so as to increase the performance at the edges of the range.

    Takes 4 arguments:

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

    random_inputs     = np.empty( (number_droplets * number_evaluations, 6), dtype=np.float32 )
    random_inputs[::number_evaluations, :] = scale_droplet_parameters( np.reshape( np.random.uniform( -1, 1, number_droplets*6 ),
                                                                       (number_droplets, 6) ).astype( "float32" ) )

    # Fix the particle temperature based on air temperature.
    random_inputs[::number_evaluations, 1] = (random_inputs[::number_evaluations, 3] +
                                              np.random.uniform( -3, 3, number_droplets ))

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
            random_inputs[droplet_index, :] = scale_droplet_parameters(
                np.random.uniform( -1, 1, 6 ).astype( "float32" ) )

            # Fix the particle temperature based on air temperature.
            random_inputs[droplet_index, 1] = (random_inputs[droplet_index, 3] +
                                               np.random.uniform( -3, 3 ))

    return random_inputs, random_outputs, integration_times

def create_training_file( file_name, number_droplets, user_batch_size=None, quiet_flag=True ):
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

    Takes 4 arguments:

      file_name       - Path to the file to write the droplet parameters.  This
                        is overwritten if it exists.
      number_droplets - The number of droplets to generate.
      user_batch_size - Optional batch size specifying the number of parameters
                        to generate at once.  If omitted,
      quiet_flag      - Optional flag specifying whether execution should be quiet.
                        If omitted, defaults to True.

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
             times) = create_droplet_batch( batch_size )

            inputs_outputs = np.hstack( (inputs,
                                         times.reshape( (batch_size, 1) ),
                                         outputs) )

            # Serialize the array
            inputs_outputs.tofile( output_fp )

def detect_timeline_gaps( integration_times, maximum_gap_size=MAXIMUM_DT_GAP_SIZE ):
    """
    Empirically identifies gaps in a time line so they may be explicitly
    handled.  Takes the successive differences between time line observations
    and computes a Boolean mask array identifying where gaps are believed to
    have occurred.  Also returns an array of time step indices indicating
    where gaps occur relative to arrays that have been masked.

    NOTE: This function handles identifying gaps in a time line of integration
          times, which implicitly corresponds to a sequence of evaluations, as
          gap handling is only reliable when operating on evaluations and their
          associated integration times.  Gap handling with observations requires
          directly working with de-gapped integration times which defeats the
          purpose of using the simulation times that naturally occurs with
          observations.

    Takes 2 arguments:

      integration_times - NumPy array, with N elements, containing the
                          successive differences between time line entries.
      maximum_gap_size  - Optional scalar specifying the maximum allowed size
                          between successive two time line entries.  Differences
                          larger than this value are flagged as a gap.  If
                          omitted, defaults to MAXIMUM_DT_GAP_SIZE.

    Returns 2 values:

      gap_mask    - NumPy Boolean array, with N elements, representing
                    differences that are considered to be gaps, where
                    N is the length of integration_times.

      gap_indices - NumPy array, with M elements (with M <= N),
                    specifying the indices where gaps occur.  These
                    indices correspond to
                    integration_times[~evaluations_gap_mask].

                    NOTE: This array may be shorter than gap_mask's True entries
                          count when multiple gaps occur back-to-back and
                          aggregate to fewer, larger gaps.

    """

    #
    # NOTE: Detection is only provided on integration times (read:
    #       evaluation-based gaps) to:
    #
    #         1. Be simple
    #         2. Avoid unnecessary footguns
    #
    #       Correctly handling gaps in sequences of observations, associated
    #       with simulations, is hard.  It is also requires more care when
    #       applying the gaps mask to construct an un-gapped timeline than for
    #       evaluations.  Since there is no real gain for supporting both and
    #       one is harder than the other, we make the obvious choice and only
    #       support evaluations-based gap detection.  This is motivated by the
    #       fact that it is *very* difficult to track down bugs when seemingly
    #       two identical operations (one observations-based, one
    #       evaluations-based) fail, only to realize they only differ by how the
    #       gaps were detected and applied.
    #
    #       For future archaeologists, the decision to only support integration
    #       times was realized after learning the following painful lessons:
    #
    #         - Construction of an observations-based gap mask is conceptually
    #           different to reason about, even though it is relatively easy to
    #           create from an evaluations-based mask.  Extending the first
    #           entry of the mask, and remove the last gap of each consecutive
    #           run (while eliminating gaps of length one), is all that is
    #           needed.
    #
    #           However, one has to keep in mind what is being masked as
    #           observations are used twice in each evaluation - once as an
    #           input and once as an output.  The mask is only applied to the
    #           observations as an input, unless there are back-to-back gaps
    #           in which case its output use is also masked.  This asymmetry
    #           can be difficult to keep track of during its use and is an
    #           easy source of mistakes.
    #
    #         - Because the observation-based gap mask is different, it also
    #           necessitates different gap indices as observations have one
    #           extra entry compared to evaluations.  As a result,
    #           observation-based each gap index is offset by its 1-based
    #           position (i.e. 1:N, for N gaps), relative to its
    #           evaluation-based gap index.  This change of index further
    #           complicates reasoning when troubleshooting differences between
    #           evaluations and observations.
    #
    #         - Most importantly, observation-based gap masks must be propagated
    #           throughout code that works on observation's simulation time
    #           and not integration times.  Unlike evaluation-based masks,
    #           one cannot simply mask gap-aware arrays and process them with
    #           their corresponding integration times *unless* they explicitly
    #           provide masked integration times (thereby defeating the point
    #           of using simulation time).  Without this, the calculated
    #           integration times will include the gaps.
    #
    #       Given that gaps are the result of computational resource constraints
    #       and are always shifting, it is likely that the same sharp edges will
    #       never be encountered twice unless testing is adopted.  The cases to
    #       worry about don't always occur (single gaps, runs of gaps, gaps at
    #       the beginning or end of a timeline that are ignored, etc) in data
    #       used for testing and requiring it to work in two similar, but subtly
    #       different ways is a recipe for frustration.
    #
    #       TL;DR: Friends don't let friends use observation-based gaps.
    #

    # Handle the case where we got zero integration times or fewer than two
    # times (the difference of on time is an empty array).
    if len( integration_times ) == 0:
        return np.array( [], dtype=np.bool_ ), np.array( [], dtype=np.int32 )

    # The gap threshold is set at 5 times the median so we're robust to
    # actual outliers in the observations, clamped at the maximum integration
    # time seen in NTLP simulations.  Some simulations have larger time
    # steps to begin with and are ratcheted down to small steps once a
    # sufficient number of particles are present.  Median will step down to
    # a smaller time step for particles with many observations while
    # selecting a reasonable step size for small observation counts.
    #
    # We have to be careful to handle the extreme corner cases where taking
    # the median breaks down, which is when we only have two integration times
    # from three observations.  This leaves us to handle the following three
    # cases:
    #
    #   1. Neither pair of observations has a gap
    #   2. One pair of observation has a gap
    #   3. Both pairs of observations have a gap
    #
    # Case #1 is handled normally by construction, though case #2 requires us to
    # compute the minimum integration time instead of the median, as this
    # results in the average of the two integration times and yields a large
    # threshold that doesn't detect the gap.  Taking the minimum at least flags
    # the gap.
    #
    # Note that Case #3 cannot be detected from a single pair of observations
    # alone unless both are beyond the maximum NTLP integration time.  As a
    # result, we opt for simplicity instead of requiring the caller to do
    # additional calculations to provide a better estimate (e.g. a running
    # median across all timelines).
    #
    # NOTE: This has not been tested on a wide range of simulations and
    #       may not be sufficiently robust!
    #
    if len( integration_times ) < 3:
        selection_func = np.min
    else:
        selection_func = np.median

    dt_threshold = min( 5 * selection_func( integration_times ),
                        MAXIMUM_DT_GAP_SIZE )

    # Compute where gaps are in the evaluations.  Derive gap indices
    # from the mask to handle runs of back-to-back gaps.
    evaluations_gap_mask = (integration_times > dt_threshold)
    gap_indices          = _evaluations_gap_mask_to_indices( evaluations_gap_mask )

    return evaluations_gap_mask, gap_indices

def _evaluations_gap_mask_to_indices( evaluations_gap_mask ):
    """
    Constructs the indices for inserting gap markers (e.g. NaN) into an
    evaluation array to produce a gapped time series.  This is useful for
    plotting where NaNs representing gaps result in visualizing multiple line
    segments instead of a single segment that spans multiple gaps.

    The indices returned are suitable for use with NumPy's .where() method.

    Gaps spanning multiple elements are compressed into a single gap.  Gaps
    at the beginning or end of the supplied mask are ignored as these
    are typically not useful.  As a result, using the indices calculated
    to construct a gapped time series does not result in the original
    mask and will have a different length.

    Takes 1 values:

      evaluations_gap_mask - NumPy array containing a Boolean mask specifying where
                             gaps occur within a time series.

    Return 1 value:

      gap_indices - NumPy array containing indices where gap markers
                    should be inserted to create the gapped array.

    """

    # We don't want to produce indices for gaps at the beginning or end
    # of the mask.  Find the subset of the mask to operate on that ignores
    # gaps at the edges.  Masks comprised entirely of gaps are properly
    # handled.
    mask_start_index = 0
    mask_end_index   = len( evaluations_gap_mask ) - 1

    while mask_start_index <= mask_end_index and evaluations_gap_mask[mask_start_index]:
        mask_start_index += 1
    while mask_end_index >= mask_start_index and evaluations_gap_mask[mask_end_index]:
        mask_end_index -= 1

    # Find gap runs where we have consecutive gaps in the mask.  We track the
    # index of the beginning of each run.
    #
    # NOTE: We only operate on the interior of the mask, after gaps starting
    #       at the beginning and up to a gap at the end.
    #
    gap_run_starts = []
    in_gap         = False
    for evaluation_index, is_gap in enumerate( evaluations_gap_mask[mask_start_index:(mask_end_index+1)] ):
        if is_gap and not in_gap:
            # Start of a gap.
            gap_run_starts.append( evaluation_index + mask_start_index )
            in_gap = True
        elif not is_gap and in_gap:
            # End of a gap.
            in_gap = False

    # Calculate the insertion position based on where we identified gaps.
    gap_indices = []
    for gap_index in gap_run_starts:
        # Count how many non-gap elements come before this gap run.
        number_nongaps_before = np.sum( ~evaluations_gap_mask[mask_start_index:gap_index] )
        gap_indices.append( number_nongaps_before )

    gap_indices = np.array( gap_indices, dtype=np.int32 )

    return gap_indices

def get_evaluation_column_names( evaluation_tag ):
    """
    Returns the names of evaluation column names for accessing the radii and
    temperature evaluations in a particles DataFrame.

    Note that these columns are only valid for DataFrames loaded with the same
    evaluations.  See read_particles_data() for details.

    Takes 1 argument:

      evaluation_tag - String specifying the evaluation file to get the
                       corresponding DataFrame columns.  May be specified as
                       a sequence of evaluation tags if multiple evaluations
                       are requested.

    Returns 1 value:

      column_names - A list of two DataFrame column names, the first component
                     is for radii and the second for the temperatures.  Will be
                     a sequence of pairs if evaluation_tag is a sequence.

    """

    # Column name templates.
    #
    # NOTE: This should match the style of the input BE column names.
    #
    RADII_TEMPLATE        = "output {:s} radii"
    TEMPERATURES_TEMPLATE = "output {:s} temperatures"

    # Figure out whether we need to wrap/unwrap a single pair or not.
    scalar_flag = isinstance( evaluation_tag, str )

    if scalar_flag:
        evaluation_tag = [evaluation_tag]

    # Build the column names for each tag.
    column_names = []
    for evaluation_tag in evaluation_tag:
        column_names.append( [RADII_TEMPLATE.format( evaluation_tag ),
                              TEMPERATURES_TEMPLATE.format( evaluation_tag )] )

    # Return the same data type as received from the caller.
    if scalar_flag:
        return column_names[0]
    else:
        return column_names

def get_evaluation_file_path( particle_path, evaluation_extension ):
    """
    Returns the path to the evaluation file associated with the provided
    raw particle file path.  The type of evaluation file is given by the
    file extension provided.

    Takes 2 arguments:

      particle_path        - Path to the reference, raw particle file.
      evaluation_extension - File name extension of the evaluation file
                             of interest.  May be specified as a sequence
                             of extensions if multiple evaluation files
                             are requested.

    Returns 1 value:

      evaluation_path - Path to the evaluation file.  Will be a sequence of
                        evaluation file paths if evaluation_extension is.

    """

    # Figure out whether we need to wrap/unwrap a single path or not.
    scalar_flag = isinstance( evaluation_extension, str )

    if scalar_flag:
        evaluation_extension = [evaluation_extension]

    # Add the file extensions to the particle file's path prefix.
    particle_root    = os.path.splitext( particle_path )[0]
    evaluation_paths = list( map( lambda extension: particle_root + "." + extension,
                                  evaluation_extension ) )

    # Return the same data type as received from the caller.
    if scalar_flag:
        return evaluation_paths[0]
    else:
        return evaluation_paths

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

def get_particles_index_path( particles_root ):
    """
    Gets a particle directory's index file path.

    Takes 1 argument:

      particles_root - Path to the top-level of a particle's directory.

    Returns 1 value:

      particle_index_path - Path to the particles directory's index.

    """

    return "{:s}/particles.index".format( particles_root )

def get_particles_parameter_extrema_path( particles_root ):
    """
    Gets a particle directory's parameters extrema file path.

    Takes 1 argument:

      particles_root - Path to the top-level of a particle's directory.

    Returns 1 value:

      particle_parameter_extrema_path - Path to the particles directory's
                                        parameter extrema file.

    """

    return "{:s}/particles.extrema".format( particles_root )

def get_particles_timeline_path( particles_root ):
    """
    Gets a particle directory's simulation timeline path.

    Takes 1 argument:

      particles_root - Path to the top-level of a particle's directory.

    Returns 1 value:

      particle_timeline_path - Path to the particles directory's simulation
                               timeline file.

    """

    return "{:s}/particles.timeline".format( particles_root )

def insert_timeseries_gaps( x, gap_indices, gap_value=np.nan ):
    """
    Takes an array representing a time series and inserts one or more values
    at each of the locations specified.  This is useful for inserting NaNs
    into a time series that will be plotted so lines are segmented at the
    gap location instead of connected.

    Takes 3 arguments:

      x           - NumPy array containing the input time series.  x may be
                    multi-dimensional though gaps will always be inserted along
                    the leading (first) dimension.
      gap_indices - NumPy array containing the indices in x indicating the right
                    side of the gap is located (i.e. first point after the gap).
                    May be empty if no gaps are present.  Indices larger than the
                    length of x are ignored so a cropped version of x can have
                    gaps inserted without the caller directly working with the
                    gaps.
      gap_value   - Optional scalar to represent the gap when inserted into
                    the time series.  If omitted, defaults to np.nan.

    Returns 1 value:

      gapped_x - NumPy array, with leading dimension of length N +
                 len( gap_indices ), containing a copy of x with
                 len( gap_indices )-many gap_values inserted into it at each
                 gap location.

    """

    if len( gap_indices ) == 0:
        return x

    # Only insert gaps in the portion of x we were provided.
    cropped_gap_indices = gap_indices[gap_indices < len( x )]

    # We need a floating point-valued array to insert NaNs.  Figure out
    # which data type to convert to.  This is a no-op for floating point
    # data types.
    if np.issubdtype( x.dtype, np.floating ):
        target_dtype = x.dtype
    elif x.itemsize > 4:
        target_dtype = np.float64
    else:
        target_dtype = np.float32

    #
    # NOTE: We insert along the first axis so we handle both 1D and 2D arrays.
    #
    return np.insert( x.astype( target_dtype ), cropped_gap_indices, np.nan, axis=0 )

def _read_particles_data_wrapper( particles_root, particle_ids, dirs_per_level, kwargs ):
    """
    Wrapper to read_particles_data() that can be called with a
    multiprocessing.pool().  This proxies a required keyword arguments
    dictionary, possibly empty, into the optional keyword arguments dictionary
    that the wrapped function expects.

    Takes 4 arguments:

      particles_root - Path to the top-level directory to read raw particle
                       files from.
      particle_ids   - Sequence of particle identifiers to read the associated
                       raw particle files.
      dirs_per_level - Number of subdirectories per level in the particle files
                       directory hierarchy.
      kwargs         - Dictionary of keyword arguments to supply to read_particles_data().

                       NOTE: This is required and not optional!

    Returns 1 value:

      particle_df - DataFrame with one row per particle identifier in particle_ids.
                    See read_particles_data() for details on the columns contained.

    """

    return read_particles_data( particles_root,
                                particle_ids,
                                dirs_per_level,
                                **kwargs )

def read_particles_data( particles_root, particle_ids, dirs_per_level, quiet_flag=True, cold_threshold=-np.inf, filter_be_failures=False, time_range=[-np.inf, np.inf], evaluations={}, debug_time_ranges=False ):
    """
    Reads one or more raw particle files and creates a particle-centric
    DataFrame with one row per particle read.  Each row file contains all of the
    observations stored as 1D NumPy arrays making it easy to operate on the
    entirety of individual particle's lifetimes.

    Filtering out outlier observations due to NTLP's backward Euler (BE)
    implementation is possible, discarding either or both "cold" particles and
    those that explicitly failed to converge during BE.

    NOTE: Filtering observations does not remove them from evaluations!
          To ensure that evaluations are aligned to observations, evaluations
          are computed on *unfiltered* observations and then filters out
          outliers.  To properly remove them from analysis, these need to be
          pruned from the raw particle files which is currently not implemented.

          As a result, this is provided as a way to simplify downstream analysis
          by factoring out outlier removal.

    Takes 9 arguments:

      particles_root     - Path to the top-level directory to read raw particle
                           files from.
      particle_ids       - Sequence of particle identifiers to read the associated
                           raw particle files.
      dirs_per_level     - Number of subdirectories per level in the particle files
                           directory hierarchy.
      quiet_flag         - Optional flag specifying whether DataFrame creation
                           should be quiet or not.  If omitted, defaults to True and
                           creation does not generate outputs unless an error
                           occurs.
      cold_threshold     - Optional floating point value specifying the lower bound
                           for particle temperatures before they're considered "cold"
                           and should be trimmed from the DataFrame.  Cold particles
                           are invalid from the NTLP simulation's standpoint but were
                           not flagged as such.  If omitted, all particles are retained.
      filter_be_failures - Optional boolean flag specifying whether observations
                           resulting from a backward Euler (BE) failure should be trimmed
                           from the DataFrame.  BE failures result in an output that
                           matches its input.  If omitted, all particles are retained.
      time_range         - Optional sequence of two values specifying the start and
                           end simulation times to crop particle observations to.
                           Observations outside of the window are discarded.  If
                           omitted, defaults to [-np.inf, np.inf] and all observations
                           are retained.
      evaluations        - Optional dictionary containing a map of tags to evaluation
                           file extensions to read and incorporate into the DataFrame.
                           If omitted, defaults to an empty dictionary and no
                           evaluations are processed.

                           NOTE: If an evaluation file is missing, DataFrame
                                 construction is aborted so the issue can be
                                 remedied.

      debug_time_ranges  - Optional Boolean flag specifying whether the cropping
                           observations to time_range should be print diagnostic
                           information to standard out or not.  If omitted, it
                           defaults to False and no information is displayed.

    Returns 1 value:

      particle_df - DataFrame with one row per particle identifier in particle_ids.
                    Contains the following columns:

                      "number evaluations":     Number of evaluations
                                                (input/output pairs for the
                                                particle.  This is the length of
                                                each of the 1D array columns.
                      "birth time":             Simulation time when this
                                                particle was created.
                      "death time":             Simulation time when this
                                                particle was destroyed.
                      "times":                  1D NumPy array of Simulation
                                                times for each of the
                                                observations.
                      "integration times":      1D NumPy array of the timestep
                                                size.
                      "number gaps":            Number of gaps in the
                                                evaluations.  This specifies the
                                                number of elements in the "gap
                                                indices" column.
                      "gap indices":            Indices specifying where a gap
                                                occurred.  These are zero-based
                                                and indicate the first
                                                evaluation after each gap.
                      "number be failures":     Non-negative number of backward
                                                Euler failures across a
                                                particle's observations.  See
                                                "be statuses" column for what
                                                constitutes a failure.
                      "be statuses":            1D NumPy array of integral
                                                backward Euler (BE) status
                                                codes, one per observation.
                                                Status is encoded as one or
                                                more bit flags:

                                                  1 - BE was skipped due to
                                                      the particle having a
                                                      negative qinf value
                                                  2 - The particle had an
                                                      equilibrium radius going
                                                      into BE
                                                  4 - BE's Gauss-Newton solve
                                                      failed
                                                  8 - Levenberg-Marquardt
                                                      failed to converge
                                                 16 - BE produced an invalid
                                                      value for either the
                                                      radius or temperature.
                                                      This represents the
                                                      cases where 1) NaNs for
                                                      either radius or
                                                      temperature, 2) radius or
                                                      temperature are negative,
                                                      or 3) radius is too large
                                                      (e.g. 10mm or larger).
                                                 32 - The particle was located
                                                      in the bottom 25% of the
                                                      vertical domain (1st
                                                      quartile)
                                                 64 - The particle was located
                                                      in the 25-50% of the
                                                      vertical domain (2nd
                                                      quartile)
                                                128 - The particle was located
                                                      in the 50-75% of the
                                                      vertical domain (3rd
                                                      quartile)
                                                256 - The particle was located
                                                      in the top 25% of the
                                                      vertical domain (4th
                                                      quartile)

                                                BE failed when the status
                                                has any of bits 1 (negative
                                                qinf), 4 (Levenberg-Marquardt
                                                failed), or 5 (invalid output)
                                                set.

                      "input be radii":         1D NumPy array of particle
                                                radii before backward Euler, in
                                                meters.
                      "output be radii":        1D NumPy array of particle
                                                radii after backward Euler, in
                                                meters.
                      "input be temperatures":  1D NumPy array of particle
                                                temperatures before backward
                                                Euler, in Kelvin.
                      "output be temperatures": 1D NumPy array of particle
                                                temperatures after backward
                                                Euler, in Kelvin.
                      "salt solutes":           1D NumPy array of the particle's
                                                salt solute, in kilograms.
                      "air temperatures":       1D NumPy array of air
                                                temperatures at the particle's
                                                positions, in Kelvin.
                      "relative humidities":    1D NumPy array of fractional
                                                relative humidities, in the
                                                range of [0, 1.5].
                      "air densities":          1D NumPy array of air density
                                                at the particle's positions,
                                                in kilogram/meter^3.

                    If evaluations is supplied then two additional columns
                    per evaluation tag (key) are present.  These tags
                    are substituted into the following templates to create
                    the additional column names:

                     "output <tag> radii":        1D NumPy array of particle
                                                  radii after <tag>'s evaluation,
                                                  in meters.
                     "output <tag> temperatures": 1D NumPy array of particle
                                                  temperatures after <tag>'s
                                                  evaluation, in Kelvin.

    """

    #
    # NOTE: We jump through hoops to get Pandas to respect our wishes regarding
    #       column data types throughout this function.  In particular, we have
    #       columns with dtype=object so we can store 1D NumPy arrays.  For the
    #       most part this works fine, though we have to work around Panda's
    #       "help" when we encounter arrays with a single element as Pandas
    #       will promote the scalar into a Series object of a single value,
    #       ignoring the parent DataFrame's data type.  While this functions
    #       *mostly* like a 1D NumPy array, it is subtly and very much not
    #       the same in that a single DataFrame row's column is now represented
    #       as a scalar meaning it's shape becomes zero-length and would require
    #       special treatment in *all* downstream callers.
    #
    #       To work around this "help", we must create a 2D array, shaped 1x1,
    #       so that Pandas doesn't promote it to a Series object (and, thus,
    #       creating a Series within a Series) and then abuse the fact that we
    #       can reshape a NumPy array by manipulating its .shape attribute.
    #       This results in a 1D array with a single element, though with much
    #       more effort than simply assigning the result of an array
    #       expression...
    #

    # Guard against an evaluations dictionary that contains the backward Euler
    # tag.  We remove this tag since there isn't a distinct evaluations file
    # for it and allows the caller to reuse the same dictionary during reading
    # and analysis.
    local_evaluations = evaluations.copy()
    if BE_TAG_NAME in local_evaluations.keys():
        del local_evaluations[BE_TAG_NAME]

    # Guard against duplicate particle identifiers that cause duplicate indices
    # in the DataFrame and breaks the logic throughout this function that
    # assumes every particle is unique.  Uniqify the identifiers while
    # maintaining the original sort order, in case it is important.
    particle_ids                            = np.array( particle_ids )
    unique_particle_ids, identifier_indices = np.unique( particle_ids, return_index=True )
    identifier_indices                      = np.sort( identifier_indices )
    particle_ids                            = particle_ids[identifier_indices]

    # Names of the observation columns.  Each of these stores an array of
    # observations that could be trimmed due to cold particles, if requested.
    #
    # NOTE: These do not yet include the evaluation outputs.  Those columns are
    #       added later.
    #
    BE_RADII_NAME, BE_TEMPERATURES_NAME = get_evaluation_column_names( BE_TAG_NAME )
    OBSERVATION_NAMES = [
        "times",
        "integration times",
        "be statuses",
        BE_RADII_NAME,
        BE_TEMPERATURES_NAME,
        "input be temperatures",
        "input be radii",
        "salt solutes",
        "air temperatures",
        "relative humidities",
        "air densities"
    ]

    # Rearrange the evaluations map into arrays of tags and extensions for
    # easier access.
    evaluation_tags       = list( local_evaluations.keys() )
    evaluation_extensions = list( local_evaluations.values() )

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

    # Create the columns we'll populate below with the correct data types.
    # Scalar columns must specify their data type (implicitly by initializing it
    # with a NumPy array) otherwise they'll default to np.float64.
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
    particles_df["number evaluations"]    = zeros_int32
    particles_df["birth time"]            = zeros_float32
    particles_df["death time"]            = zeros_float32
    particles_df["times"]                 = pd.Series( dtype=object )
    particles_df["integration times"]     = pd.Series( dtype=object )
    particles_df["number gaps"]           = zeros_int32
    particles_df["gap indices"]           = pd.Series( dtype=object )
    particles_df["number be failures"]    = zeros_int32
    particles_df["be statuses"]           = pd.Series( dtype=object )
    particles_df["input be radii"]        = pd.Series( dtype=object )
    particles_df[BE_RADII_NAME]           = pd.Series( dtype=object )
    particles_df["input be temperatures"] = pd.Series( dtype=object )
    particles_df[BE_TEMPERATURES_NAME]    = pd.Series( dtype=object )
    particles_df["salt solutes"]          = pd.Series( dtype=object )
    particles_df["air temperatures"]      = pd.Series( dtype=object )
    particles_df["relative humidities"]   = pd.Series( dtype=object )
    particles_df["air densities"]         = pd.Series( dtype=object )

    # Add the evaluation columns.  Each of them are an array like the other
    # observation columns.
    #
    # NOTE: We add these last as they're not (necessarily) more important than
    #       the original observations.
    #
    evaluation_column_names = get_evaluation_column_names( evaluation_tags )
    for column_names_pair in evaluation_column_names:
        evaluation_radii_name, evaluation_temperatures_name = column_names_pair

        particles_df[evaluation_radii_name]        = pd.Series( dtype=object )
        particles_df[evaluation_temperatures_name] = pd.Series( dtype=object )

        OBSERVATION_NAMES.extend( [evaluation_radii_name, evaluation_temperatures_name] )

    # Start processing particles.
    for particle_id in particle_ids:
        particle_path = get_particle_file_path( particles_root, particle_id, dirs_per_level )

        # Get each of the observations.
        observations_fp32  = _read_raw_particle_data( particle_path )
        observations_int32 = observations_fp32.view( np.int32 )

        if not quiet_flag:
            print( "Processing {:d}".format( int( particle_id ) ) )

        # Build this particle's timeline so we can keep the observations of
        # interest.
        #
        # NOTE: This mask includes all observation's whose inputs occur within
        #       the range, admitting the last observation's corresponding output
        #       to fall just outside of the range.
        #
        observations_mask_all = ((observations_fp32[:, ParticleRecord.TIME_INDEX] >= time_range[0]) &
                                 (observations_fp32[:, ParticleRecord.TIME_INDEX] <= time_range[1]))

        # Setup our indexing to pull the correct observations for constructing
        # the particle's evaluations, as well as when we read additional
        # evaluations.  Applying Boolean mask arrays is a heavyweight operation
        # relative to slicing (3-10x slower, depending on the number of
        # observations) so we only use them when we have to.
        #
        # NOTE: Supporting a timeline window makes things complicated as we have
        #       to both slice the observations and, separately, mask them.
        #       We cannot compose both of these operations together into a single
        #       slice-like access as applying a Boolean mask will shrink the
        #
        if (observations_mask_all.sum() == observations_fp32.shape[0]):
            # Slices encode indexing operations (i.e.. ":", ":-1", and "1:").
            observations_mask_all     = slice( None )
            observations_mask_outputs = slice( 1, None )
            observations_mask_inputs  = slice( None, -1 )
            evaluations_mask_all      = slice( None )

            # We have one evaluation for every output observation.
            number_evaluations = observations_fp32.shape[0] - 1
        else:
            # We're cropping available observations to a window of time, so
            # we must use Boolean mask arrays that are each slightly altered
            # to represent both a subset in time and a subset of that subset
            # (e.g ignoring the first or last cropped value).
            #
            # The inputs and the outputs masks are strict subsets of the all
            # mask, each retaining one fewer observation than the all mask.  The
            # input mask is aligned at the beginning of the all mask, while the
            # output mask is shifted forward by one observation and is aligned
            # at the end of the all mask.
            #
            # Evaluations are always time aligned with the available
            # observations' timeline and *NOT* the currently requested time
            # window.  This means we simply ignore the first observation.
            observations_mask_inputs  = observations_mask_all.copy()
            observations_mask_outputs = observations_mask_all.copy()
            evaluations_mask_all      = observations_mask_all[:-1].copy()

            # Determine if any of the observations fall within the time range
            # requested.  This particle's observations could lie entirely
            # out of the range which means we already have the masks of
            # interest that are all False.
            observations_included_mask = np.where( observations_mask_all )[0]
            if len( observations_included_mask ) > 0:
                # Ranges that do not include the final observation need to
                # be extended by one observation so as to include it as the
                # original mask's last input observation's output.  This
                # covers the crop cases of (oo, B] and [A, B].
                if observations_included_mask[-1] != (len( observations_mask_all ) - 1):
                    observations_mask_all[observations_included_mask[-1] + 1]     = True
                    observations_mask_outputs[observations_included_mask[-1] + 1] = True
                else:
                    # We're handling the crop case of [A, oo) where a subset
                    # of observations ending with the last available has been
                    # requested.  Remove the final included observation from
                    # the input mask.
                    observations_mask_inputs[observations_included_mask[-1]] = False

                # Regardless of the crop case ((oo, B], [A, oo), or [A, B]) we
                # need to remove the first included observation from the output
                # mask.  This "shifts" the input observations to the outputs.
                observations_mask_outputs[observations_included_mask[0]] = False

                # Print the masks if we're debugging.  Note that evaluations
                # match the output observations so we align their mask so its
                # easier to see.
                if debug_time_ranges:
                    print( "\n"
                           "All:         {}\n",
                           "Inputs:      {}\n",
                           "Outputs:     {}\n",
                           "Evaluations:       {}\n".format(
                               observations_mask_all,
                               observations_mask_inputs,
                               observations_mask_outputs,
                               evaluations_mask_all ) )

            # We have one evaluation for every output observation.
            number_evaluations = max( 0, observations_mask_outputs.sum() )

        # We have one evaluation for every output observation.
        particles_df.at[particle_id, "number evaluations"] = number_evaluations

        # Keep track of when this particle existed regardless of whether we
        # retain any of it's observations.
        particles_df.at[particle_id, "birth time"] = observations_fp32[0,  ParticleRecord.TIME_INDEX]
        particles_df.at[particle_id, "death time"] = observations_fp32[-1, ParticleRecord.TIME_INDEX]

        # Use the particle's observation times to compute the evaluations'
        # integration times.
        particles_df.at[particle_id, "integration times"] = np.diff( observations_fp32[observations_mask_all, ParticleRecord.TIME_INDEX] )

        # We assume that we have all observations for this particle.
        particles_df.at[particle_id, "number gaps"] = 0
        particles_df.at[particle_id, "gap indices"] = np.array( [] )

        # Skip particles that don't have any evaluations.  These are either
        # because they have zero observations or because they have one
        # evaluation with a large dt that implies there was a gap between the
        # first and last observation.  In both cases we retain the birth and
        # death times and wipe out the individual observations.
        #
        # NOTE: Observe the lack of subscript on the "integration times" as
        #       Pandas has converted the singleton array to a scalar!
        #
        if ((particles_df.at[particle_id, "number evaluations"] == 0) or
            (particles_df.at[particle_id, "number evaluations"] == 1 and
             particles_df.at[particle_id, "integration times"] >= MAXIMUM_DT_GAP_SIZE)):

            particles_df.at[particle_id, "number evaluations"] = 0
            particles_df.at[particle_id, "number be failures"] = 0

            for observation_name in OBSERVATION_NAMES:
                particles_df.at[particle_id, observation_name] = np.array( [] )

            continue

        # Backward Euler failures for this particle.
        particles_df.at[particle_id, "be statuses"]           = observations_int32[observations_mask_inputs, ParticleRecord.BE_STATUS_INDEX]
        particles_df.at[particle_id, "number be failures"]    = BEStatus.get_failure_mask( particles_df.at[particle_id, "be statuses"] ).sum()

        # Observed radii.
        particles_df.at[particle_id, "input be radii"]        = observations_fp32[observations_mask_inputs,  ParticleRecord.RADIUS_INDEX]
        particles_df.at[particle_id, BE_RADII_NAME]           = observations_fp32[observations_mask_outputs, ParticleRecord.RADIUS_INDEX]

        # Observed temperatures.
        particles_df.at[particle_id, "input be temperatures"] = observations_fp32[observations_mask_inputs,  ParticleRecord.TEMPERATURE_INDEX]
        particles_df.at[particle_id, BE_TEMPERATURES_NAME]    = observations_fp32[observations_mask_outputs, ParticleRecord.TEMPERATURE_INDEX]

        # Background parameters.
        particles_df.at[particle_id, "salt solutes"]          = observations_fp32[observations_mask_inputs, ParticleRecord.SALT_SOLUTE_INDEX]
        particles_df.at[particle_id, "air temperatures"]      = observations_fp32[observations_mask_inputs, ParticleRecord.AIR_TEMPERATURE_INDEX]
        particles_df.at[particle_id, "relative humidities"]   = observations_fp32[observations_mask_inputs, ParticleRecord.RELATIVE_HUMIDITY_INDEX]
        particles_df.at[particle_id, "air densities"]         = observations_fp32[observations_mask_inputs, ParticleRecord.AIR_DENSITY_INDEX]

        # The particle's time line.
        particles_df.at[particle_id, "times"]                 = observations_fp32[observations_mask_inputs, ParticleRecord.TIME_INDEX]

        # Get the evaluation file paths to read and the corresponding columns.
        evaluation_file_paths = get_evaluation_file_path( particle_path,
                                                          evaluation_extensions )
        for evaluation_file_path, column_names_pair in zip( evaluation_file_paths,
                                                            evaluation_column_names ):

            (evaluation_radii_name,
             evaluation_temperatures_name) = column_names_pair
            (evaluation_radii,
             evaluation_temperatures)      = _read_particle_evaluation_data( evaluation_file_path )

            # Keep the evaluations of interest.
            particles_df.at[particle_id, evaluation_radii_name]        = evaluation_radii[evaluations_mask_all]
            particles_df.at[particle_id, evaluation_temperatures_name] = evaluation_temperatures[evaluations_mask_all]

        # Work around Pandas being "helpful" by converting a single element
        # array into a scalar.  Without this, two observation/single evaluation
        # particles explode on use as the assumption that we have arrays of
        # evaluations is fundamental to this DataFrame.
        if len( particles_df.at[particle_id, "times"].shape ) == 0:
            for observation_name in OBSERVATION_NAMES:
                particles_df.at[particle_id, observation_name] = np.array( [[particles_df.at[particle_id, observation_name]]] )
                particles_df.at[particle_id, observation_name].shape = [1]

        # Drop the evaluations that have crazy integration times under the
        # assumption that we don't have a complete trace of this particle (some
        # observations occurred on a MPI rank whose outputs were not processed).
        # This manifests as outputs that appear to have a huge integration time
        # because the intermediate observations are missing.
        #
        # NOTE: We use number of masked evaluations to determine whether we
        #       update this particle's arrays, and not the length of the
        #       indices, as the former accounts for gaps at the beginning/end of
        #       the time line that are omitted from the gap indices.
        #
        gaps_mask, gap_indices    = detect_timeline_gaps( particles_df.at[particle_id, "integration times"] )
        number_masked_evaluations = gaps_mask.sum()
        number_gap_indices        = len( gap_indices )

        if number_masked_evaluations > 0:
            # We have at least one gap that we need to remove.  Since we're have
            # evaluations with an input and an output we can simply drop those
            # with large integration times.
            #
            # NOTE: We could do this earlier if we absolutely cared about every
            #       last ounce of performance by carefully managing which
            #       observations could be used as inputs but not outputs.  It is
            #       much easier to do it this way, albeit at the expense of
            #       shifting all of this particle's arrays over an element (or
            #       few).
            #
            #       This shouldn't happen on fully evaluated datasets so we
            #       don't consider this a priority just yet.
            #

            # Track how many gaps, and where they're located, so downstream
            # users can correctly segment the ungapped regions.  We take care to
            # adjust the gap locations to account for the fact that the original
            # evaluations containing them will be removed below.
            particles_df.at[particle_id, "number gaps"] = number_gap_indices
            if number_gap_indices > 1:
                particles_df.at[particle_id, "gap indices"] = gap_indices
            elif number_gap_indices == 1:
                # Work around Pandas' "help".  Note that we don't need to adjust
                # the locations since we only have a single gap.
                particles_df.at[particle_id, "gap indices"]       = np.array( [gap_indices[0]] )
                particles_df.at[particle_id, "gap indices"].shape = [1]

            # Move the birth time forward if we have a gap between the first two
            # observations, which manifests itself in the first evaluation.
            new_start_index = np.argmin( gaps_mask )
            if new_start_index > 0:
                particles_df.at[particle_id, "birth time"] = particles_df.at[particle_id, "times"][new_start_index]

            # Move the death time backward if we have a gap between the last two
            # observations, which manifests itself in the last evaluation.
            new_end_index = -np.argmin( gaps_mask[::-1] )
            if new_end_index < 0:
                particles_df.at[particle_id, "death time"] = particles_df.at[particle_id, "times"][new_end_index]

            # Reduce the count of BE failures that occurred on an observation
            # preceding a gap.  When this occurs the output observation
            # associated with the failure is in the gap and the evaluation
            # incorrectly combines the BE failure input with the first
            # observation after the gap.
            #
            # NOTE: We do this before removing the gaps so as to correctly count
            #       how many evaluations are removed below.
            #
            number_gapped_be_failures = (gaps_mask & BEStatus.get_failure_mask( particles_df.at[particle_id, "be statuses"] )).sum()
            particles_df.at[particle_id, "number be failures"] -= number_gapped_be_failures

            # Get all of the evaluation columns as a single, flat list.
            all_evaluation_column_names = sum( evaluation_column_names, [] )

            # Remove the gaps from every observation.
            for observation_name in OBSERVATION_NAMES:
                # Avoid a "chicken and the egg" problem and skip evaluation
                # columns.  These already had gaps removed prior to their
                # creation and don't need them removed again.
                if observation_name in all_evaluation_column_names:
                    continue

                particles_df.at[particle_id, observation_name] = particles_df.at[particle_id, observation_name][~gaps_mask]

                # Once again, work around Pandas being "helpful" by converting
                # a single element array into a scalar.  Without this, indexing
                # this particle's columns will explode.
                if number_masked_evaluations == (particles_df.at[particle_id, "number evaluations"] - 1):
                    # Work around Pandas' "help".
                    particles_df.at[particle_id, observation_name] = np.array( [[particles_df.at[particle_id, observation_name]]] )
                    particles_df.at[particle_id, observation_name].shape = [1]

            # Finally, reduce the observation count by the gaps we've removed.
            particles_df.at[particle_id, "number evaluations"] -= number_masked_evaluations

            # Move to the next particle if we eliminated all of the evaluations.
            if particles_df.at[particle_id, "number evaluations"] == 0:
                continue

        # Provide a safety net against incorrect indexing.  We want to
        # immediately know about broken indexing instead of randomly tripping up
        # on a particle or, worse, never noticing and having subtly incorrect
        # analysis.
        #
        # This ensures we have matching counts across all observation arrays.
        # We take care to track the individual counts to aide in debugging
        # should we explode.
        #
        # NOTE: We have to do this after removing gaps since evaluations are
        #       created with gaps already removed.  Otherwise we're guaranteed
        #       to have mismatched array sizes for particles with gaps.
        #
        array_sizes        = set()
        observation_counts = []
        for observation_name in OBSERVATION_NAMES:
            if len( particles_df.at[particle_id, observation_name].shape ) > 0:
                # 1D array.
                observation_count = particles_df.at[particle_id, observation_name].shape[0]
            else:
                # Scalar that slipped through our conversion contortions.
                observation_count = 1

            array_sizes.add( observation_count )
            observation_counts.append( observation_count )

        if len( array_sizes ) != 1:
            raise ValueError( "Array size mismatch for particle id {:d} ({}, {})!".format(
                particle_id,
                observation_counts,
                OBSERVATION_NAMES ) )

        # Excise failed backward Euler (BE) evaluations if requested.  This
        # removes evaluations whose inputs match their BE outputs.
        if filter_be_failures and particles_df.at[particle_id, "number be failures"] > 0:

            # BE failures typically occur at the beginning of a particle's life
            # though have been observed throughout simulations.
            be_failure_mask = BEStatus.get_failure_mask( particles_df.at[particle_id, "be statuses"] )

            # Count the number of successive failures at the beginning so we can
            # adjust the particle's birth time accordingly.
            if (particles_df.at[particle_id, "be statuses"][0] & BEStatus.FAILED_UPDATE) > 0:

                # Compute the index of the first non-failed observation so we
                # can use it's simulation time as the start of this particle's
                # life.
                number_starting_failures = 1
                for observation_index in range( 1, particles_df.at[particle_id, "number evaluations"] ):
                    if (particles_df.at[particle_id, "be statuses"][observation_index] & BEStatus.FAILED_UPDATE) == 0:
                        break

                    number_starting_failures += 1

                # Check the corner case where *ALL* of the observations are BE
                # failures and excising them will result in an empty particle.
                # We retain the original birth and death times in that case so
                # as to be consistent with the time range filtering from above.
                if number_starting_failures != particles_df.at[particle_id, "number evaluations"]:
                    particles_df.at[particle_id, "birth time"] = particles_df.at[particle_id, "times"][number_starting_failures]

            # Reduce the number of observations we have.
            particles_df.at[particle_id, "number evaluations"] -= be_failure_mask.sum()
            particles_df.at[particle_id, "number be failures"] -= be_failure_mask.sum()

            # Excise the failed BE observations.
            for observation_name in OBSERVATION_NAMES:
                particles_df.at[particle_id, observation_name] = particles_df.at[particle_id, observation_name][~be_failure_mask]

        # Excise cold observations if requested.  These are an artifact of
        # backward Euler that produces extremely cold temperatures when the
        # process should have failed.
        if particles_df.at[particle_id, "input be temperatures"].min() < cold_threshold:

            # Cold particles occur once toward the beginning of their lifetime
            # where their temperature plunges for one (maybe two) timesteps.

            # Trim one additional point after the cold particle recovers.
            NUMBER_PADDING_POINTS = 1

            # Locate where the cold particle occurs.
            cold_observation_index = particles_df.at[particle_id, "input be temperatures"].argmin()
            trim_index             = cold_observation_index + 1 + NUMBER_PADDING_POINTS

            # Subtract the trimmed observations and pretend that the particle
            # sprang into existence after our cut.
            particles_df.at[particle_id, "number evaluations"] -= trim_index + 1
            particles_df.at[particle_id, "birth time"]          = particles_df.at[particle_id, "times"][trim_index]

            # Excise the observations in question.
            for observation_name in OBSERVATION_NAMES:
                particles_df.at[particle_id, observation_name] = particles_df.at[particle_id, observation_name][trim_index:]

    return particles_df

def _read_particle_evaluation_data( evaluation_file_path ):
    """
    Reads an evaluation file and returns the evaluations' radii and temperatures
    as two NumPy arrays.  Separate arrays are returned to simplify assignment of
    the resulting data into DataFrame columns.

    NOTE: It is assumed that the contents of the evaluations file is sorted by
          simulation time!  This means the arrays returned are *NOT* guaranteed
          to be time aligned with raw particle observations.

    Takes 1 argument:

      evaluation_file_path - Path to evaluation file to read.

    Returns 2 values:

      evaluation_radii        - NumPy array, sized number_time_steps x 1, containing
                                the evaluation radii.
      evaluation_temperatures - NumPy array, sized number_time_steps x 1, containing
                                the evaluation temperatures.

    """

    # Open the file and read it in big chunk (possibly with multiple filesystem reads)
    # and let NumPy create an array from the buffer.
    #
    # NOTE: We transpose the data so it's in column-major order to match the
    #       access patterns typically encountered by the caller (e.g. operating
    #       on a single field for all observations).
    #
    # NOTE: Manually reading the data and creating an array from a buffer is
    #       *significantly* faster than simply using .fromfile().  Depending on
    #       the number of observations read, and the underlying file system,
    #       this can be 3-10x faster.
    #
    with open( evaluation_file_path, "rb" ) as evaluation_fp:
        evaluation_data = evaluation_fp.read()
    #
    # NOTE: Evaluation data are already sorted by time!
    #
    evaluation_data = np.frombuffer( evaluation_data, dtype=np.float32 ).reshape( (-1, 2) ).copy( order="F" )

    # Split the data into two separate arrays since we're assigning these to
    # Pandas Series.
    #
    # NOTE: These are really views on evaluation_data from above, so
    #       the underlying storage order matters.
    #
    evaluation_radii        = evaluation_data[:, 0]
    evaluation_temperatures = evaluation_data[:, 1]

    return (evaluation_radii, evaluation_temperatures)

def _read_raw_particle_data( particle_path ):
    """
    Reads a raw particle file and returns its contents as an array of 32-bit
    floating point values, sorted by the simulation time field.  The integer
    field values are not cast as floating point values but rather are simply
    interpreted as floating point values at the bit level and can be accessed
    by applying a view to the returned array.

    Takes 1 argument:

      particle_path - Path to the raw particle file to read.

    Returns 1 value:

      observations_fp32 - NumPy array, shaped number_observations x
                          ParticleRecord.SIZE, containing the particle's
                          observations interpreted as 32-bit floating point
                          values.  Since two of the columns are actually 32-bit
                          integer values, access to those columns must be
                          performed via a .view() on the returned array.

                          See ParticleRecord for details on the array's columns.

    """

    # Open the file and read it in big chunk (possibly with multiple filesystem reads)
    # and let NumPy create an array from the buffer.
    #
    # NOTE: We transpose the data so it's in column-major order to match the
    #       access patterns typically encountered by the caller (e.g. operating
    #       on a single field for all observations).
    #
    # NOTE: Manually reading the data and creating an array from a buffer is
    #       *significantly* faster than simply using .fromfile().  Depending on
    #       the number of observations read, and the underlying file system,
    #       this can be 3-10x faster.
    #
    with open( particle_path, "rb" ) as particle_fp:
        particle_data = particle_fp.read()
    observations_fp32 = np.frombuffer( particle_data, dtype=np.float32 ).reshape( -1, ParticleRecord.SIZE ).copy( order="F" )

    # Sort the observations by time.
    #
    # NOTE: This is needed in case this particle moved across MPI ranks during
    #       the simulation and was handled by more than one rank.  In that case
    #       the observations within each rank are sorted but we could get
    #       windows of observations out of order depending on how the resulting
    #       trace files were post-processed into raw particle files.
    #
    sorted_indices    = np.argsort( observations_fp32[:, ParticleRecord.TIME_INDEX] )
    observations_fp32 = observations_fp32[sorted_indices, :]

    return observations_fp32

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

def write_particles_index( particles_index_path, particle_ids, merge_flag=False ):
    """
    Creates or updates a particles index file using the supplied particle
    identifiers.  The unique set of particle identifiers provided is written to
    the supplied particles index.  If requested, an existing index file's
    contents can be merged with the supplied identifiers facilitating
    incremental updates to a particles directory hierarchy.

    This function handles concurrent access by locking the particle index.

    Takes 3 arguments:

      particles_index_path - Path to the particles index file to write.  This
                             file is overwritten when not merging particle
                             indices.
      particle_ids         - Array of particle identifiers to write to the
                             index.  Can also be a list of arrays.
      merge_flag           - Optional Boolean specifying whether particle_ids
                             should be merged into an existing index (True) or a
                             new index created from particle_ids (False).  If
                             omitted, defaults to False and a new index is
                             created, overwriting particles_index_path if it
                             exists.

    Returns 1 value:

      unique_particle_ids - A NumPy array of unique particle identifiers in the
                            index.  If merge_flag is True then this also
                            includes the identifiers that were previously in
                            particles_index_path.

    """

    # Help the caller and convert a list of arrays into a single array.
    if isinstance( particle_ids, list ):
        particle_ids = np.hstack( particle_ids )

    # Open the file for read/write, creating it if necessary, but not truncating
    # it if it exists.  This allows us to lock the file once while we update the
    # index which is more efficient when we are merging indices with the
    # contents of the index (i.e. two lock/unlocks).
    with open( particles_index_path, "ab+" ) as index_fp:

        try:
            fcntl.flock( index_fp.fileno(), fcntl.LOCK_EX )

            # Seek to the beginning and read the index if we're merging
            # our particles in.
            if merge_flag:

                # Determine the file's length.  If it's non-zero then
                # we have an existing set of particles to merge with.
                index_fp.seek( 0, os.SEEK_END )
                if index_fp.tell() > 0:
                    index_fp.seek( 0 )

                    # Read the current index.
                    #
                    # NOTE: We have at least one particle (for well-formed
                    #       indices) since the file size is non-zero.
                    #
                    existing_particle_ids = np.fromfile( index_fp,
                                                         dtype=np.int32 )

                    # Add our particles to the end of the current index.
                    particle_ids = np.hstack( [existing_particle_ids,
                                               particle_ids] )

            # De-duplicate the particles.
            unique_particle_ids = np.unique( particle_ids )

            # Overwrite the index with the unique particles.
            index_fp.seek( 0 )
            index_fp.truncate()
            unique_particle_ids.tofile( index_fp )

            # Attempt to flush any buffered data to disk before releasing our
            # lock.
            index_fp.flush()

        except IOError as e:
            raise RuntimeError( "Failed to lock '{:s}': {:s}".format(
                particles_index_path,
                str( e ) ) )

        finally:
            # Ensure the file is unlocked for the next read/writer.
            fcntl.flock( index_fp.fileno(), fcntl.LOCK_UN )

    return unique_particle_ids
