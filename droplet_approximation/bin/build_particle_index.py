#!/usr/bin/env python3

"""
Processes a NTLP trace file and extracts metadata to create a particle index, the
observations' parameter extrema, and an estimate of the simulation timeline.
This is intended to be run once per NTLP trace and write to a new particle
directory tree so that additional post-processing can be launched to populate
the tree.

Care is taken to minimize the calculations performed so that this script runs
as fast as I/O allows, as there is often a lot of data to process.

The particle index is key for downstream processing and is required to exist for
some scripts to even be queued for future launch (e.g. creating evaluations from
the observations).  Parameter extrema are useful for understanding the bounds of
the observations to inform future development and analysis decisions.  Creating
a simulation timeline simplifies analysis as post-processed observations are
often too large to hold in memory, preventing creation of the timeline from the
data loaded during analysis.

Usage: build_particle_index.py <trace_path> <particle_root>

"""

import fcntl
import getopt
import os
import sys
from typing import Optional

import numpy as np

import droplet_approximation
from droplet_approximation import ParticleRecord

# Read data in 16 MiB chunks so we can let the system perform efficient reads.
READ_CHUNK_SIZE = 16 * 1024 * 1024

def merge_array_with_lock( output_path: str, values: np.ndarray, extrema_flag: Optional[bool] =False ):
    """
    Merges the contents of a file with the array supplied and overwrites the
    file with the unique, sorted combination.

    Takes 3 arguments:

        output_path  - Path to the file to update.  It is overwritten if it
                       exists, or created otherwise.
        values       - NumPy array of values to write.
        extrema_flag - Optional Boolean specifying that values is 2D, shaped N x 2,
                       and contains parameter extrema.  If omitted, defaults
                       to False and values must be a 1D array.

    Returns nothing.

    """

    bytes_per_element = values.dtype.itemsize

    # Open the file for appending so we can read and write it.  This creates it
    # if it doesn't exist.
    with open( output_path, "ab+" ) as output_fp:
        # Acquire an exclusive lock during the entire operation.
        fcntl.flock( output_fp.fileno(), fcntl.LOCK_EX )

        try:
            # Read the entirety of file into memory.
            output_fp.seek( 0 )
            file_data = output_fp.read()

            existing_values = np.array( [], dtype=values.dtype )

            if len( file_data ) > 0:
                # Ensure the file has an integral number of elements in it.
                number_trailing_bytes = len( file_data ) % bytes_per_element
                if number_trailing_bytes != 0:
                    print( "'{:s}' is not a multiple of {:d} bytes ({:d} % {:d} = {:d})!".format(
                        output_path,
                        len( file_data ),
                        len( file_data ),
                        bytes_per_element,
                        number_trailing_bytes ) )

                    # Truncate the partial bytes to force an integral number of
                    # elements.
                    file_data = file_data[:len(file_data) - number_trailing_bytes]

                if file_data:
                    existing_values = np.frombuffer( file_data, dtype=values.dtype )

            # Combine the old and new values depending on whether we have
            # 2D extrema or a 1D array of values.
            if extrema_flag:
                values          = values.reshape( (-1, 2) )
                existing_values = existing_values.reshape( (-1, 2) )

                # Deal with the case where we're creating a new extrema file.
                if existing_values.shape[0] == 0:
                    combined_values = values
                else:
                    combined_values       = np.empty_like( values )
                    combined_values[:, 0] = np.minimum( values[:, 0],
                                                        existing_values[:, 0] )
                    combined_values[:, 1] = np.maximum( values[:, 1],
                                                        existing_values[:, 1] )
            else:
                combined_values = np.unique( np.concatenate( [existing_values,
                                                              values] ) )

            # Wipe out the existing file contents, if any.
            output_fp.seek( 0 )
            output_fp.truncate()

            # Write the combined values to disk.
            combined_values.flatten().tofile( output_fp )

        finally:
            # Ensure that the lock is released regardless of how we leave this
            # block.
            fcntl.flock( output_fp.fileno(), fcntl.LOCK_UN )

def process_observations( int32_records: np.ndarray, fp32_records: np.ndarray ):
    """
    Processes a chunk of particle observations to compute:

      1. The unique particle identifiers
      2. The parameter extrema
      3. The simulation times

    NOTE: This does not exactly compute the extrema for the integration times
          parameter as that requires a more computationally expensive approach
          that dramatically slows this function without realistically improving
          the reported integration times.  See the note in this routine for more
          details.

    Takes 2 arguments:

      int32_records - NumPy array, shaped N x 2, containing an 32-bit integer
                      view of the particle observations to process.  This is
                      aliased to the same memory as fp32_records.
      fp32_records  - NumPy array, shaped N x 2, containing an 32-bit floating
                      point view of the particle observations to process.  This
                      is aliased to the same memory as int32_records.

    Returns 3 values:

      particle_ids      - NumPy array containing the unique, sorted, 32-bit
                          integer particle identifiers in int32_records.
      parameter_extrema - NumPy array, shaped ParticleRecord.SIZE x 2, containing
                          the minimum and maximums for each of ParticleRecord's
                          fields.  The PARTICLE_ID and BE_STATUS_INDEX entries
                          are unused and are returned with undefined values.
      simulation_times  - NumPy array containing the unique, sorted, 32-bit
                          floating point simulation times in fp32_records.

    """

    # Compute the unique values.
    particle_ids = np.unique( int32_records[:, 0] )

    # Recover this chunk's simulation timeline.
    simulation_times = np.unique( fp32_records[:, ParticleRecord.TIME_INDEX] )

    # Track the extrema of each non-temporal parameter.
    parameter_extrema = np.empty( (ParticleRecord.SIZE, 2), dtype=np.float32 )

    parameter_extrema[ParticleRecord.RADIUS_INDEX, :]            = (fp32_records[:, ParticleRecord.RADIUS_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.RADIUS_INDEX].max())
    parameter_extrema[ParticleRecord.TEMPERATURE_INDEX, :]       = (fp32_records[:, ParticleRecord.TEMPERATURE_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.TEMPERATURE_INDEX].max())
    parameter_extrema[ParticleRecord.SALT_SOLUTE_INDEX, :]       = (fp32_records[:, ParticleRecord.SALT_SOLUTE_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.SALT_SOLUTE_INDEX].max())
    parameter_extrema[ParticleRecord.AIR_TEMPERATURE_INDEX, :]   = (fp32_records[:, ParticleRecord.AIR_TEMPERATURE_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.AIR_TEMPERATURE_INDEX].max())
    parameter_extrema[ParticleRecord.RELATIVE_HUMIDITY_INDEX, :] = (fp32_records[:, ParticleRecord.RELATIVE_HUMIDITY_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.RELATIVE_HUMIDITY_INDEX].max())
    parameter_extrema[ParticleRecord.AIR_DENSITY_INDEX, :]       = (fp32_records[:, ParticleRecord.AIR_DENSITY_INDEX].min(),
                                                                    fp32_records[:, ParticleRecord.AIR_DENSITY_INDEX].max())

    # Calculating extrema for integration time is tricky since we need two
    # consecutive observations from the same particle.  Looping through each
    # particle available in the chunk gets us all but one potential integration
    # time - the one that uses the last observation from each particle in the
    # previous chunk and computes the difference from the first observation of
    # the corresponding particle in this chunk.
    #
    # Since we aren't exactly computing the extrema for the integration times we
    # also approximate the integration times so as to greatly accelerate this
    # computation.  Looping through every particle's integration times and
    # computing the extrema is quite time consuming and ends up being the
    # dominant computation for the entire script (by about 100x).  Instead we
    # simply take the difference of the unique times which should result in the
    # same extrema.  This approach is acceptable as speed is significantly more
    # important and the extrema reported here are typically padded before use,
    # so being off by 1e-6 is acceptable.
    #
    dts = np.diff( simulation_times )
    parameter_extrema[ParticleRecord.TIME_INDEX, :] = (dts.min(), dts.max())

    return particle_ids, parameter_extrema, simulation_times

def process_trace_file( trace_path: str, particles_index_path: str, particles_extrema_path: str,
                        particles_timeline_path: str, chunk_size: int ):
    """
    Processes a NTLP trace file, extracts key metadata, and writes it to disk
    for additional post-processing.  Particle observations are scanned for the
    following:

      1. Unique list of particle identifiers to construct an index
      2. Parameter extrema calculated across all observations
      3. Unique simulation times to build a simulation time line

    Processing is done in chunks, with the size specified by the caller, so as
    to balance I/O efficiency against memory footprint.  Chunk sizes are rounded
    down to the nearest particle observation record size.

    Takes 5 arguments:

      trace_path              -
      particles_index_path    -
      particles_extrema_path  -
      particles_timeline_path -
      chunk_size              -

    Returns 1 value

    Args: XXX

        trace_path (str): Path to input binary file
        particles_index_path (str): Path to output binary file
        particles_extrema_path (str): Path to output binary file
        chunk_size (int): Size of chunks to read at once (in bytes)
        fields_per_record (int): Number of 32-bit fields per record

    """

    unique_particle_ids     = np.array( [], np.int32 )
    unique_simulation_times = np.array( [], np.float32 )

    record_size = ParticleRecord.SIZE_BYTES

    # Ensure chunk_size is a multiple of record_size for clean reading
    chunk_size = (chunk_size // record_size) * record_size
    if chunk_size == 0:
        chunk_size = record_size

    # Initialize the extrema.
    #
    # NOTE: This includes extrema for non-parameters (e.g. particle identifier)
    #       that we don't actually compute.  This is done so we can use
    #       ParticleRecord's indices without complicating the code.
    #
    parameter_extrema = np.empty( (ParticleRecord.SIZE, 2), dtype=np.float32 )
    parameter_extrema[:, 0] =  np.inf
    parameter_extrema[:, 1] = -np.inf

    with open( trace_path, "rb" ) as trace_fp:
        print( "Processing {:s}".format( trace_path ) )

        records_processed = 0
        chunks_processed  = 0

        while True:
            # Read chunk of observations.
            chunk = trace_fp.read( chunk_size )
            if not chunk:
                break

            # Only process complete records.  XXX: .read() should only return full records.
            complete_bytes = (len( chunk ) // record_size) * record_size
            if complete_bytes > 0:
                # Create 32-bit floating point view of the data read and
                # then create a transposed copy so we have the records in
                # column-major order.  This provides unit stride access when
                # operating across all records.
                fp32_view    = np.frombuffer( chunk[:complete_bytes], dtype='<f4' )
                fp32_records = fp32_view.reshape( -1, ParticleRecord.SIZE ).copy( order="F" )

                # Create 32-bit integer view of the records.
                #
                # NOTE: These are also in column-major order.
                #
                int32_records = fp32_records.view( dtype=np.int32 )

                # Parse the observations.
                (particle_ids,
                 current_parameter_extrema,
                 simulation_times) = process_observations( int32_records,
                                                           fp32_records )

                # Update our extrema.
                np.minimum( parameter_extrema[:, 0],
                            current_parameter_extrema[:, 0],
                            out=parameter_extrema[:, 0] )
                np.maximum( parameter_extrema[:, 1],
                            current_parameter_extrema[:, 1],
                            out=parameter_extrema[:, 1] )

                # Update our particles list and timeline with this chunk's
                # unique values.
                unique_particle_ids     = np.unique( np.concatenate( (unique_particle_ids,
                                                                      particle_ids) ) )
                unique_simulation_times = np.unique( np.concatenate( (unique_simulation_times,
                                                                      simulation_times) ) )

                records_processed += len( int32_records )
                chunks_processed  += 1

            # Handle incomplete record at end of chunk
            remaining_bytes = len( chunk ) - complete_bytes
            if remaining_bytes > 0:
                # Seek back to start of incomplete record
                trace_fp.seek( -remaining_bytes, 1 )

            # Progress indicator
            if chunks_processed % 1 == 0:
                print( "Processed {:d} records, found {:d} unique values".format(
                    records_processed,
                    len( unique_particle_ids )
                ) )

    print( "\n"
           "Total records processed: {:d}\n"
           "Unique values found:     {:d}".format(
               records_processed,
               len( unique_particle_ids ) )  )

    # Write out our contributions to the index, the parameter extrema, and
    # timeline files.
    #
    # NOTE: The extrema array is sized to simplify indexing but we only care
    #       about the particle characteristics and not the placeholders for
    #       particle identifiers and backward Euler status.
    #
    merge_array_with_lock( particles_index_path, unique_particle_ids )
    merge_array_with_lock( particles_extrema_path,
                           parameter_extrema[ParticleRecord.TIME_INDEX:, :],
                           extrema_flag=True )
    merge_array_with_lock( particles_timeline_path, unique_simulation_times )

    return True

def main( argv ):
    """
    Main entry point for the script.  Parses the command line and launches
    process_trace_file().

    Takes 1 argument:

      argv - Command line components, including the script's name as the first
             component.

    Returns 1 value:

      exit_code - Integral code to exit the script with.

    """

    usage_string = \
        "Usage: {:s} <trace_path> <particles_root>".format(
            argv[0] )

    # Parse any options (none currently used, but following the pattern)
    try:
        options, arguments = getopt.getopt( argv[1:], "" )
    except getopt.GetoptError as error:
        print( "{:s}\n"
               "\n"
               "{:s}".format(
                   error,
                   usage_string ) )
        return 1

    for option, option_value in options:
        raise NotImplementedError( "Unhandled option '{:s}'!".format( option ) )

    NUMBER_REQUIRED_ARGUMENTS = 2

    number_arguments = len( arguments )

    if number_arguments != NUMBER_REQUIRED_ARGUMENTS:
        print( usage_string )
        return 1

    trace_path     = arguments[0]
    particles_root = arguments[1]

    print( "Trace:     {:s}".format( trace_path ) )
    print( "Particles: {:s}".format( particles_root ) )

    # Create the particles root if needed.
    os.makedirs( particles_root, exist_ok=True )

    # Get the paths to the files we're creating/updating.
    particles_index_path    = droplet_approximation.get_particles_index_path( particles_root )
    particles_extrema_path  = droplet_approximation.get_particles_parameter_extrema_path( particles_root )
    particles_timeline_path = droplet_approximation.get_particles_timeline_path( particles_root )

    # Process the file.
    success_flag = process_trace_file(
        trace_path,
        particles_index_path,
        particles_extrema_path,
        particles_timeline_path,
        READ_CHUNK_SIZE
    )

    if not success_flag:
        return 1

    return 1

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
