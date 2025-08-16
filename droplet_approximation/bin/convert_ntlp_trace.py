#!/usr/bin/env python3

"""
Trace file processing script that distributes work across multiple processes.

This script processes trace files containing fixed-size records using multiprocessing.
Each record consists of 9 32-bit float values (36 bytes total). The script can be
run as part of a larger distributed job system where multiple instances process
different subsets of the total data.

Usage:
    script.py [-t] <particles_root> <number_jobs> <job_index> <file1> [<file2> [...]]

    -t: Test mode - print configuration without executing

"""

import multiprocessing
import os
import sys
from typing import List, Tuple

import numpy as np

import droplet_approximation
import droplet_approximation.data

# Convenience variables to (slightly) reduce the name verbosity.  These specify
# size of each NTLP observation both in values and bytes.
VALUES_PER_RECORD = droplet_approximation.data.ParticleRecord.SIZE.value
RECORD_SIZE       = droplet_approximation.data.ParticleRecord.SIZE_BYTES.value

# Read a large number of records at once to minimize read overhead.
BATCH_SIZE = 1024*1024

# Number of directories per level in the particles hierarchy.  A branching
# factor of 256 provides roughly 65K directories to store particles in so as to
# not store too many particles in any one.  This is tuned for hierarchies
# containing millions of particles.
DIRS_PER_LEVEL = 256

def main( arguments: List[str] ):
    """
    Main entry point for the trace processing script.

    Parses command line arguments, calculates work distribution, and either
    prints test information or launches worker processes.

    Takes 1 argument:

      arguments - List of strings provided on the command line to parse.  This
                  includes any options as well as the required arguments.

    Returns nothing.

    """

    # Parse the provided command line arguments.
    (test_mode,
     particles_root,
     number_jobs,
     job_index,
     trace_files) = parse_arguments( arguments )

    # Calculate the total number of records across all input files.
    total_records = calculate_total_records( trace_files )

    # Determine records this job instance should process.
    (record_start_index,
     record_end_index) = calculate_job_range( total_records, number_jobs, job_index )

    # Calculate how work will be distributed among worker processes.
    number_processes = 1 #multiprocessing.cpu_count()
    worker_ranges    = divide_range_among_workers( record_start_index,
                                                   record_end_index,
                                                   number_processes )

    # Show this job's configuration.
    print_summary( trace_files,
                   number_jobs,
                   job_index,
                   record_start_index,
                   record_end_index,
                   number_processes,
                   worker_ranges,
                   particles_root )

    # Nothing left to do if we're only testing.
    if test_mode:
        return

    # Create a pool of worker processes to handle the work.
    if number_processes == 1:
        particle_ids_list = worker_process( 0,
                                            worker_ranges[0][0],
                                            worker_ranges[0][1],
                                            trace_files,
                                            particles_root )
    else:
        with multiprocessing.Pool( processes=number_processes ) as pool:

            # Create the arguments for each of the workers.
            args_list = []
            for worker_id, (start_record_index, end_record_index) in enumerate( worker_ranges ):
                args_list.append( (worker_id,
                                   start_record_index,
                                   end_record_index,
                                   trace_files,
                                   particles_root) )

            # Do particle extraction.
            particle_ids_list = pool.starmap( worker_process, args_list )

            pool.close()
            pool.join()

    # Merge the processed particle identifiers into the index, creating
    # it if necessary.
    particles_index = droplet_approximation.get_particles_index_path( particles_root )
    droplet_approximation.write_particles_index( particles_index,
                                                 particle_ids_list,
                                                 merge_flag=True )

def parse_arguments( arguments: List[str] ) -> Tuple[bool, str, int, int, List[str]]:
    """
    Parse and validate command line arguments.  Exits if provided with an
    incorrect number of arguments or given invalid arguments.

    Raises ValueError if arguments are supplied of the wrong type and cannot be
    converted to the expected type.

    Takes 1 argument:

      arguments - Array of command line arguments to parse.  This does not
                  include the executing script's name.

    Returns 5 values:

      testing_flag   - Boolean indicating if test mode was requested.
      particles_root - Path to top-level of a particles directory.
      number_jobs    - Total number of job instances.
      job_index      - Index of this job instance.
      trace_files    - List of input trace file paths.

    """

    # Validate minimum number of required arguments.
    if len( arguments ) < 4:
        print( "Error: Insufficient arguments.\n"
               "\n"
               "Usage: script.py [-t] <particles_root> <number_jobs> <job_index> <file1> [<file2> [...]]",
               file=sys.stderr )
        sys.exit( 1 )

    # Check for optional test mode flag.
    testing_flag = False
    if arguments[0] == '-t':
        testing_flag = True

        # Remove the flag from remaining arguments.
        arguments = arguments[1:]

    # Extract and validate arguments.
    particles_root = arguments[0]
    number_jobs    = int( arguments[1] )
    job_index      = int( arguments[2] )
    trace_files    = arguments[3:]

    # Validate job parameters.
    if number_jobs < 1:
        print( "Number of jobs must be positive ({:d})!".format(
            number_jobs ),
               file=sys.stderr )
        sys.exit( 1 )

    if job_index < 0 or job_index >= number_jobs:
        print( "The job index ({:d}) must be between 0 and {:d}!".format(
            job_index,
            number_jobs - 1 ),
              file=sys.stderr )
        sys.exit( 1 )

    # Verify all trace files exist.
    for file_path in trace_files:
        if not os.path.isfile( file_path ):
            print( "Trace file '{:s}' not found!".format( file_path ),
                   file=sys.stderr )
            sys.exit( 1 )

    return testing_flag, particles_root, number_jobs, job_index, trace_files

def get_file_size( file_path: str ) -> int:
    """
    Get the size of a file in bytes.

    Takes 1 argument:

      file_path - Path to the file.

    Returns 1 value:

      file_size - Size of the file in bytes.

    """

    return os.path.getsize( file_path )

def human_readable_count( count: int ) -> str:
    """
    Converts a count into a human-readable string.  For example, 1250000000
    becomes "1.25 billion".

    Takes 1 argument:

      count - Count to transform into a string.

    Returns 1 value:

      count_str - Human-readable string approximating count.

    """

    if count > 10**9:
        count_str = "{:.2f} billion".format( count / 10**9 )
    elif count > 10**6:
        count_str = "{:.2f} million".format( count / 10**6 )
    elif count > 10**3:
        count_str = "{:.2f} thousand".format( count / 10**3 )
    else:
        count_str = str( count )

    return count_str

def human_readable_size( number_bytes: int ) -> str:
    """
    Converts a byte size into a human-readable string.  For example,
    1250000000 becomes "1.16 GiB".

    Takes 1 argument:

      number_bytes - Number of bytes to transform into a string.

    Returns 1 value:

      size_str - Human-readable string approximating number_bytes.

    """

    if number_bytes > 1024**4:
        size_str = "{:.2f} TiB".format( number_bytes / 1024**4 )
    elif number_bytes > 1024**3:
        size_str = "{:.2f} GiB".format( number_bytes / 1024**3 )
    elif number_bytes > 1024**2:
        size_str = "{:.2f} MiB".format( number_bytes / 1024**2 )
    elif number_bytes > 1024**1:
        size_str = "{:.2f} KiB".format( number_bytes / 1024**1 )
    else:
        size_str = "{:d} bytes".format( number_bytes )

    return size_str

def calculate_total_records( trace_files: List[str] ) -> int:
    """
    Calculate the total number of records across all trace files.

    Raises ValueError if any of the supplied files has a partial record.

    Takes 1 argument:

      trace_files - List of trace file paths.

    Returns 1 value:

      number_records - Total number of records across all files.

    """

    total_file_sizes = 0

    for trace_file in trace_files:
        file_size = get_file_size( trace_file )

        if file_size % RECORD_SIZE > 0:
            raise ValueError( "'{:s}' has a partial trailing record ({:d} bytes)!".format(
                trace_file,
                file_size % RECORD_SIZE ) )

        total_file_sizes += file_size

    number_records = total_file_sizes // RECORD_SIZE

    return number_records

def calculate_job_range( number_total_records: int, number_jobs: int, job_index: int ) -> Tuple[int, int]:
    """
    Calculate the range of records this job instance should process.

    Takes 3 arguments:

      number_total_records - Total number of records across all files.
      number_jobs          - Total number of job instances.
      job_index            - Index of this job instance.

    Returns 2 values:

      start_record_index - Global record index to start at.
      end_record_index   - Global record index to end with, exclusively.

    """

    records_per_job   = number_total_records // number_jobs
    records_remainder = number_total_records % number_jobs

    # Compute the start and stop indices for this job.  Excess records are given
    # to the first few processes.
    if job_index < records_remainder:
        start_record_index = job_index * (records_per_job + 1)
        end_record_index   = start_record_index + records_per_job + 1
    else:
        start_record_index = job_index * records_per_job + records_remainder
        end_record_index   = start_record_index + records_per_job

    return start_record_index, end_record_index

def divide_range_among_workers( start_record_index: int, end_record_index: int,
                                number_workers: int ) -> List[Tuple[int, int]]:
    """
    Divide a range of records among worker processes.

    Takes 3 arguments:

      start_record_index - Starting record index for this job
      end_record_index   - Ending record index for this job (exclusive)
      number_workers     - Number of worker processes

    Returns 1 value:

      worker_ranges - Sequence of two element tuples containing the global start
                      and end record indices for each worker.

    """

    number_records     = end_record_index - start_record_index
    records_per_worker = number_records // number_workers
    remainder          = number_records % number_workers

    current_start_index = start_record_index

    # Build a list of (start, end) pairs, one per worker.
    worker_ranges = []
    for worker_index in range( number_workers ):

        # The first N workers get an extra record if the number of workers
        # doesn't evenly divide the records, where N is the remainder.
        number_worker_records = records_per_worker + int( worker_index < remainder )
        worker_end_index      = current_start_index + number_worker_records

        # Only add workers that have records to process.  This handles the case
        # where we have more workers than records.
        if number_worker_records > 0:
            worker_ranges.append( (current_start_index, worker_end_index) )

        current_start_index = worker_end_index

    return worker_ranges

def record_index_to_file_position( global_record_index: int, trace_files: List[str] ) -> Tuple[int, int]:
    """
    Convert a global record index to file index and byte offset.

    Raises ValueError if the requested record index is beyond the last record in
    the final trace file.

    Takes 2 arguments:

      global_record_index - Global record index across all files.
      trace_files         - List of trace file paths.

    Returns 2 values:

      file_index  - Index into trace_files where global_record_index is located.
      byte_offset - Byte offset into trace_files[file_index] where
                    global_record_index is located.

    """

    # Walk through the trace files and identify which contains the record of
    # interest.  Keep track of the current number of records seen in all of
    # the previous files.
    number_cumulative_records = 0
    for file_index, file_path in enumerate( trace_files ):

        file_size           = get_file_size( file_path )
        number_file_records = file_size // RECORD_SIZE

        # Is the requested record in this file?
        if global_record_index < (number_cumulative_records + number_file_records):
            record_file_offset = global_record_index - number_cumulative_records
            byte_offset        = record_file_offset * RECORD_SIZE

            return file_index, byte_offset

        number_cumulative_records += number_file_records

    raise ValueError( "Requested record index exceeds total records available ({:d} > {:d}.".format(
        global_record_index,
        number_cumulative_records ) )

def print_summary( trace_files: List[str], number_jobs: int,
                   job_index: int, job_start_record_index: int, job_end_record_index: int,
                   number_processes: int, worker_ranges: List[Tuple[int, int]],
                   particles_root: str ):
    """
    Print a comprehensive summary of the processing configuration for test mode.

    Takes 8 arguments:

      trace_files            - List of particle trace file paths to extract
                               particle observations from.
      number_jobs            - Total number of job instances.
      job_index              - Index of this job instance.
      job_start_record_index - Starting global record index for this job.
      job_end_record_index   - Ending global record index for this job.
      number_processes       - Number of worker processes.
      worker_ranges          - List of global record range pairs to process, one
                               per worker.  These record indices will be a
                               subset of [job_start_record_index,
                               job_end_record_index].
      particles_root         - Path to top-level of a particles directory.

    Returns nothing.

    """

    # Map out the files.  Build a map from file names to a tuple containing
    # their size in bytes and records.
    trace_file_map       = {}
    for trace_path in trace_files:
        trace_size_bytes     = get_file_size( trace_path )
        number_trace_records = trace_size_bytes // RECORD_SIZE

        trace_file_map[trace_path] = (trace_size_bytes, number_trace_records)

    # Compute the number of records across all files.
    total_number_records = sum( map( lambda pair: pair[1], trace_file_map.values() ) )

    # Print the overview of the conversion.
    print( "Configuration Summary:\n"
           "\n"
           "  Output directory:      {:s}\n"
           "  Traces:\n"
           "\n"
           "    Number files:        {:d}\n"
           "    Number observations: {:s} ({:d})\n"
           "    Size:                {:s} ({:d} bytes)\n"
           "\n"
           "  {:d} workers:\n"
           "\n"
           "    Batch size:          {:d} observations\n"
           "    Batches/worker:      {:d} ({:d} observations)\n"
           "    Data read/worker:    {:s} ({:d} bytes)"
           "\n".format(
               particles_root,
               len( trace_files ),
               human_readable_count( total_number_records ),
               total_number_records,
               human_readable_size( total_number_records * RECORD_SIZE ),
               total_number_records * RECORD_SIZE,
               number_processes,
               BATCH_SIZE,
               total_number_records // number_processes // BATCH_SIZE,
               total_number_records // number_processes,
               human_readable_size( total_number_records // number_processes * RECORD_SIZE ),
               total_number_records // number_processes * RECORD_SIZE ) )

    # Print details on the trace files being processed.
    print( "Trace files map:"
          "\n")

    for trace_number, trace_path in enumerate( trace_files, 1 ):
        trace_size_bytes    = get_file_size( trace_path )
        number_trace_records = trace_size_bytes // RECORD_SIZE

        print( "  {:d}. {:s} ({:d} bytes, {:d} observation{:s})".format(
            trace_number,
            trace_path,
            trace_size_bytes,
            number_trace_records,
            "" if number_trace_records == 1 else "s" ) )

    print()

    # Print details on what each worker will process.
    print( "Job #{:d} of {:d} processing {:d} record{:s} ({:d} - {:d}) has {:d} worker assignment{:s}:"
           "\n".format(
               job_index + 1,
               number_jobs,
               job_end_record_index - job_start_record_index,
               "" if (job_end_record_index - job_start_record_index) == 1 else "s",
               job_start_record_index,
               job_end_record_index,
               len( worker_ranges ),
               "" if len( worker_ranges ) == 1 else "s" ) )

    for worker_number, (worker_start_record_index, worker_end_record_index) in enumerate( worker_ranges, 1 ):
        number_records = worker_end_record_index - worker_start_record_index
        number_batches = (number_records + BATCH_SIZE - 1) // BATCH_SIZE

        print( "  {:d}. Record indices {:d} - {:d} ({:d} total in {:d} batch{:s}) from:"
               "\n".format(
            worker_number,
            worker_start_record_index,
            worker_end_record_index,
            number_records,
            number_batches,
            "" if number_batches == 1 else "es" ) )

        # Show file span for this worker.
        start_file_index, start_byte_offset = record_index_to_file_position( worker_start_record_index,
                                                                             trace_files )
        end_file_index, end_byte_offset     = record_index_to_file_position( worker_end_record_index - 1,
                                                                             trace_files )

        # Add one record to the end so we display a record range, rather than
        # the offsets of the first and last records read.  [start_index,
        # end_index] is more useful for verification and review than "reading at
        # indices start_index to end_index" which has an implicit off-by-one.
        end_byte_offset += RECORD_SIZE

        # Are we working within a single file?
        if start_file_index == end_file_index:
            print( "    {:s}: Bytes [{:d} - {:d}]".format(
                trace_files[start_file_index],
                start_byte_offset,
                end_byte_offset ) )
        else:
            # This work spans two or more files.
            for processed_file_index in range( start_file_index,
                                               end_file_index + 1 ):
                file_path       = trace_files[processed_file_index]
                file_size_bytes = get_file_size( file_path )

                # We have three cases to consider:
                #
                #   1. This is the first file and we're processing to the end of
                #      it.
                #   2. This is the last file and we're processing from the
                #      beginning of it.
                #   3. We're in between the first and last files and are
                #      processing the entirety of it.
                #
                if processed_file_index == start_file_index:
                    print( "    {:s}: Bytes [{:d} - {:d}]".format(
                        file_path,
                        start_byte_offset,
                        file_size_bytes ) )
                elif processed_file_index == end_file_index:
                    print( "    {:s}: Bytes [{:d} - {:d}]".format(
                        file_path,
                        0,
                        end_byte_offset ) )
                else:
                    print( "    {:s}: Bytes [{:d} - {:d}]".format(
                        file_path,
                        0,
                        file_size_bytes ) )

        print()

def worker_process( worker_id: int, start_record_index: int, end_record_index: int,
                    trace_files: List[str], particles_root: str ):
    """
    Main worker process function that reads particle observations from one or
    more files and writes out raw particles trajectories to the supplied
    particles directory.  The specified range of records is located in the
    provided trace files and may span multiple trace files so as to facilitate
    flexible observation processing.

    Takes 5 arguments:

      worker_id          - Unique identifier for this worker process.
      start_record_index - Starting global record index for this worker
      end_record_index   - Ending global record index for this worker (XXX: exclusive?).
      trace_files        - List of trace file paths.
      particles_root     - Path to top-level of a particles directory.

    Returns 1 value:

      unique_particle_ids - NumPy array of unique particle identifiers processed.

    """

    # List of NumPy arrays containing particle identifiers seen in each batch.
    particle_ids_list = []

    current_record_index = start_record_index
    while current_record_index < end_record_index:
        # Determine how many records to process in this batch.  Take care to not
        # go beyond the last record this worker was assigned.
        batch_size = min( BATCH_SIZE,
                          end_record_index - current_record_index )

        # Read a batch of records as a NumPy array.
        particle_observations = read_record_batch( current_record_index, batch_size, trace_files )

        # Write these observations out to the raw particle file that owns them.
        particle_ids = \
            droplet_approximation.convert_NTLP_trace_array_to_particle_files( particle_observations,
                                                                              particles_root,
                                                                              DIRS_PER_LEVEL )
        particle_ids_list.append( particle_ids )

        current_record_index += batch_size

    # Return an array of unique particle identifiers.  We will have duplicates
    # when a particle's observations span multiple batches.
    return np.unique( np.hstack( particle_ids_list ) )

def read_record_batch( start_record_index: int, batch_size: int,
                       trace_files: List[str] ) -> np.ndarray:
    """
    Read a batch of records efficiently, handling file boundaries.

    Takes 3 arguments:

      start_record_index - Starting global record index for the batch.
      batch_size         - Number of records to read.
      trace_files        - List of trace file paths.

    Returns 1 value:

      batch_data - NumPy array of shape (batch_size, VALUES_PER_RECORD)
                   containing the records as 32-bit floating point values.

    """

    # Allocate our 2D output array, one row per record in this batch.  We copy
    # our data into this array so we can correctly handle reading records that
    # span across multiple files.
    batch_data = np.empty( (batch_size, VALUES_PER_RECORD), dtype=np.float32 )

    number_records_read  = 0
    current_record_index = start_record_index

    while number_records_read < batch_size:
        # Find which file and position for current record.
        file_index, byte_offset = record_index_to_file_position( current_record_index,
                                                                 trace_files )

        # Calculate how many contiguous records we can read from this file.
        file_size                 = get_file_size( trace_files[file_index] )
        records_remaining_in_file = (file_size - byte_offset) // RECORD_SIZE
        number_records_needed     = batch_size - number_records_read
        number_records_to_read    = min( records_remaining_in_file, number_records_needed )

        # Read contiguous chunk from this file.
        number_bytes_to_read = number_records_to_read * RECORD_SIZE
        file_data            = read_file_chunk( trace_files[file_index], byte_offset, number_bytes_to_read )

        # Convert bytes to NumPy array of float32 values.
        # Reshape from flat array to (records_to_read, VALUES_PER_RECORD).
        chunk_array = np.frombuffer( file_data, dtype=np.float32 ).reshape(
            number_records_to_read, VALUES_PER_RECORD )

        # Copy into batch_data array.
        batch_data[number_records_read:number_records_read + number_records_to_read] = chunk_array

        number_records_read  += number_records_to_read
        current_record_index += number_records_to_read

    return batch_data

def read_file_chunk( file_path: str, byte_offset: int, number_bytes: int ) -> bytes:
    """
    Read a number of bytes from a file starting at a given offset.

    Takes 3 arguments:

        file_path    - Path to the file to read.
        byte_offset  - Starting byte position in file_path.
        number_bytes - Number of bytes to read.

    Returns 1 value:

        bytes_read - Bytes object containing the raw data read from file_path.

    """

    with open( file_path, "rb" ) as file_handle:
        file_handle.seek( byte_offset )

        return file_handle.read( number_bytes )

if __name__ == "__main__":
    main( sys.argv[1:] )
