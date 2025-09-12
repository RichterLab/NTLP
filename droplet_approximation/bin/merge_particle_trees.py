#!/usr/bin/env python3

"""
Script to merge contents of two particle directory trees.  Raw particle
observations are read from the source tree, merged with any existing
observations from the destination tree and then written back to the destination
tree.  Evaluations in the destination tree are removed to ensure mismatched
observations and evaluations do not occur.

There are several key things to note in this implementation:

1. Newly merged observations are sorted by time to aide in debugging and
   diagnostics.  This is allowed per the raw particle "format" but can result
   in different files being copied should the destination directory not exist.

2. Evaluations are not handled as simply merging two gapped evaluation time
   series will maintain the gap restart evaluations which will be different when
   the merge eliminates a gap.  While we could simply refuse to merge iterative
   evaluations, this provides an easy way to accidentally merge the wrong thing
   and introduce an nearly impossible to find error which impacts all downstream
   analysis.

3. This script is written to support concurrent merging of multiple source
   trees into the same destination tree!  Locking is used to allow simultaneous
   merging of multiple subsets of a particle found in different NTLP traces.

"""

import fcntl
import getopt
import glob
import os
import multiprocessing
import sys
from typing import List, Tuple

import numpy as np

import droplet_approximation
from droplet_approximation import ParticleRecord

# This script aims to provide a key tool in pipelines to incrementally process
# large particle traces so as to balance the time it takes to generate data for
# analysis against deadlines as well as the underlying storage required.
#
# We intentionally implement a "copy" of just raw particle files and force
# the re-evaluation of the resulting particles for two reasons:
#
#   1. Merging iterative evaluation files that were created from gapped
#      observations can result in different time series than if iterative
#      evaluation is applied after the merge.  Should the merge remove a gap
#      that was previously present, the merged evaluation would still reflect
#      that gap rather than the ungapped evaluation.
#
#   2. As a result from #1, this provides the flexibility for batching up
#      multiple merges into a single, staging tree that is evaluated once
#      as well as iterative merging and (redundantly) evaluating.
#
# While we could whitelist the preservation of direct (non-iterative)
# evaluations, that would provide an opportunity to accidentally miss an
# evaluation to retain or, worse, accidentally include an iterative evaluation
# that is now silently incorrect and forever impacts downstream analysis.  This
# approach aims to minimize self-inflicted wounds at the relatively minimal cost
# of re-evaluating everything.

# Number of directories per level in the particles hierarchy.  A branching
# factor of 256 provides roughly 65K directories to store particles in so as to
# not store too many particles in any one.  This is tuned for hierarchies
# containing millions of particles.
DIRS_PER_LEVEL = 256

class PartialSourceObservationsError( ValueError ):
    """
    Specialized exception indicating the source particle observations contain an incomplete observation.
    """

class PartialDestinationObservationsError( ValueError ):
    """
    Specialized exception indicating the destination particle observations contain an incomplete observation.
    """

class CouldNotRemoveEvaluationsError( RuntimeError ):
    """
    Specialized exception indicating a failure to remove an particle evaluations file.
    """

class FailedToWriteObservationsError( RuntimeError ):
    """
    Specialized exception indicating a failure to write particle observations.
    """

def merge_raw_particle_file( source_particle_path: str, destination_particle_path: str ):
    """
    Merges observations from one raw particle file into another particle file.
    The source observations are combined with the destination observations,
    sorted by simulation time, and then overwrite the destination particle
    file.

    An exclusive lock on destination_particle_path is acquired before merging
    observations into it to prevent race conditions allowing multiple callers
    to safely execute concurrently.

    NOTE: This function assumes that the parent directory of
          destination_particle_path exists!

    May raise the following exceptions:

      PartialSourceObservationsError       source_particle_path contains an
                                           incomplete observation
      PartialDestinationObservationsError  destination_particle_path contains an
                                           incomplete observation
      IOError                              Failure to acquire an exclusive lock on
                                           destination_particle_path

    Takes 2 arguments:

      source_particle_path      - Path to a raw particle file whose contents
                                  to merge into destination_particle_path.
      destination_particle_path - Path to a raw particle file to overwrite
                                  with the combination of its contents and
                                  source_particle_path's.

    Returns nothing.

    """

    # Read the source file.  Bail if we didn't read an integral number
    # of observations.
    with open( source_particle_path, "rb" ) as source_fp:
        source_buffer = source_fp.read()

        if len( source_buffer ) % ParticleRecord.SIZE_BYTES != 0:
            raise PartialSourceObservationsError(
                "'{:s}' does not contain an integral number of observations ({:d} extra bytes)!".format(
                    source_particle_path,
                    len( source_buffer ) % ParticleRecord.SIZE_BYTES ) )

    # Get a NumPy array of observations from the buffer.  This is a 32-bit
    # floating point view since we only care about the simulation time field.
    source_observations = np.frombuffer( source_buffer,
                                         dtype="<f4" ).reshape( (-1, ParticleRecord.SIZE) )

    # Open the destination file for append so we can read its previous contents
    # and overwrite it with the combined observations.
    with open( destination_particle_path, "ab+" ) as destination_fp:
        try:
            # Prevent races where we write the wrong combined file and lose
            # observations another process is writing at the same time.
            fcntl.flock( destination_fp.fileno(), fcntl.LOCK_EX )

            # Read the file.
            destination_fp.seek( 0 )
            destination_buffer = destination_fp.read()

            if len( destination_buffer ) == 0:
                # We just created the file, so we have nothing to combine with.
                destination_observations = np.array( [],
                                                     dtype=np.float32 ).reshape( (0, ParticleRecord.SIZE) )
            elif len( destination_buffer ) % ParticleRecord.SIZE_BYTES != 0:
                # Bail if we detect an incomplete observation.
                raise PartialDestinationObservationsError(
                "'{:s}' does not contain an integral number of observations ({:d} extra bytes)!".format(
                    destination_particle_path,
                    len( destination_buffer ) % ParticleRecord.SIZE_BYTES ) )
            else:
                # Get a NumPy array that is shape compatible with the source
                # observations.
                destination_observations = np.frombuffer( destination_buffer,
                                                          dtype="<f4" ).reshape( (-1, ParticleRecord.SIZE) )

            # Concatenate the two sets of observations and sort them by time.
            all_observations           = np.concatenate( [source_observations,
                                                          destination_observations],
                                                         axis=0 )
            sorted_observation_indices = np.argsort( all_observations[:, ParticleRecord.TIME_INDEX] )
            all_observations           = all_observations[sorted_observation_indices, :]

            # Overwrite the destination file.
            destination_fp.seek( 0 )
            destination_fp.truncate()
            all_observations.tofile( destination_fp )

        except OSError as e:
            raise IOError( "Failed to lock '{:s}' - {:s}".format(
                destination_particle_path,
                str( e ) ) )

        finally:
            # Always unlock the file regardless of how the merge went, otherwise
            # we screw over others attempting to safely access this file.
            fcntl.flock( destination_fp.fileno(), fcntl.LOCK_UN )

def remove_evaluation_files( particle_path: str ) -> None:
    """
    Removes any evaluation file associated with the specified raw particle file.

    Raises CouldNotRemoveEvaluationsError if removing an evaluation file
    results in a failure.

    Takes 1 argument:

      particle_path - Path to the raw particle file whose evaluations should be
                      removed.

    Returns nothing.

    """

    # Break the raw particle file's path down into its constituent components.
    particle_directory   = os.path.dirname( particle_path )
    particle_base_name   = os.path.basename( particle_path )
    (particle_base_name,
     particle_extension) = os.path.splitext( particle_base_name )

    # Get paths of all files that have the same base name as the provided
    # raw particle file.
    evaluation_files_pattern = os.path.join( particle_directory,
                                             "{:s}.{:s}".format(
                                                 particle_base_name, "*" ) )
    matching_paths = glob.glob( evaluation_files_pattern )

    # Filter out the original file so we don't nuke it.
    evaluation_file_paths = filter( lambda path: os.path.splitext( path )[1] != particle_extension,
                                    matching_paths )

    try:
        for evaluation_file_path in evaluation_file_paths:
            os.remove( evaluation_file_path )
    except (PermissionError, FileNotFoundError) as e:
        raise CouldNotRemoveEvaluationsError( "Failed to remove {:s}: {:s}".format( evaluation_file_path, e ) )

def merge_particle_batch( particle_ids: np.ndarray, source_tree: str,
                          destination_tree: str, keep_evaluations_flag: bool ) -> Tuple[List[str], List[str]]:
    """
    Merges the specified particles' observations found in one particle
    tree with the observations in another, overwriting the destination
    observations with their combination.  Returns lists of raw particle
    files that could not be read or written.

    Takes 4 arguments:

      particle_ids          - 1D NumPy array of particle identifiers to read from
                              source_tree and merge into destination_tree.
      source_tree           - Particle tree to read particles from.
      destination_tree      - Particle tree to merge into.  Particles from
                              source_tree are merged with the particles in
                              this tree.
      keep_evaluations_flag - Boolean flag indicating whether evaluations should
                              be kept or deleted in destination_tree.

    Returns 2 values:

      read_errors  - List of pairs, particle paths and the exception's message,
                     that could not be read from source_tree.
      write_errors - List of pairs, particle paths and the exception's message,
                     that could not be written to destination_tree.

    """

    # Keep track of paths that could not be merged.
    read_errors  = []
    write_errors = []

    for particle_id in particle_ids:
        # Get the source and destination paths,
        source_path = droplet_approximation.get_particle_file_path( source_tree,
                                                                    particle_id,
                                                                    DIRS_PER_LEVEL )
        destination_path = droplet_approximation.get_particle_file_path( destination_tree,
                                                                         particle_id,
                                                                         DIRS_PER_LEVEL )

        # Nothing to do if the source file doesn't exist.
        if not os.path.exists( source_path ):
            read_errors.append( source_path )
            continue

        # Ensure the destination directory exists.
        #
        # NOTE: We do this here, once per particle, to minimize the load on the
        #       file system.
        #
        os.makedirs( os.path.dirname( destination_path ), exist_ok=True )

        try:
            # Copy the source observations into the destination tree.  Wipe out
            # all of the evaluations that exist in the destination afterwards
            # since they're no longer valid.
            merge_raw_particle_file( source_path, destination_path )

            if not keep_evaluations_flag:
                remove_evaluation_files( destination_path )

        # Map exceptions to read or write failures.
        except PartialSourceObservationsError as e:
            read_errors.append( (source_path, str( e )) )
            continue
        except (PartialDestinationObservationsError,
                CouldNotRemoveEvaluationsError,
                FailedToWriteObservationsError) as e:
            write_errors.append( (destination_path, str( e )) )
            continue

    return read_errors, write_errors

def merge_directory_trees( source_tree: str, destination_tree: str, number_processes: int, keep_evaluations_flag: bool ) -> None:
    """
    Merges all of the particles found in one particles tree into another
    tree using one or more workers for concurrency.  Workers are spawned as
    separate processes.

    Raises FileNotFoundError if the source tree's particle index does not exist.

    Takes 4 arguments:

      source_tree           - Path to the source particles tree.  All particles
                              in this tree's index are merged into
                              destination_tree.
      destination_tree      - Path to the destination particles tree.  All
                              particles in this tree that are also in
                              source_tree's index are overwritten with their
                              combined observations.  This directory is created
                              if it doesn't exist.
      number_processes      - Number of worker processes to merge source_tree
                              and destination_tree with.  Must be a positive
                              integer.
      keep_evaluations_flag - Boolean flag indicating whether evaluations should
                              be kept or deleted in destination_tree.

    Returns nothing.

    """

    # Bail if the source tree doesn't have an index.
    source_index_path = droplet_approximation.get_particles_index_path( source_tree )
    if not os.path.exists( source_index_path ):
        raise FileNotFoundError( "Particles index file ('{:s}') not found!".format(
            source_index_path ) )

    # Get all of particle identifiers in the source tree.
    unique_particle_ids = np.fromfile( source_index_path, dtype=np.int32 )

    # Bail if there isn't anything to merge.  We don't raise an exception since
    # we've technically done what was requested.
    if len( unique_particle_ids ) == 0:
        print( "No particles found in '{:s}'.  Nothing to do!".format( source_index_path ) )
        return

    # Split the available particles across the workers.
    particles_per_process = len( unique_particle_ids ) // number_processes

    # Build the workers' arguments lists.
    batch_arguments = []
    for process_index in range( number_processes ):
        start_index = process_index * particles_per_process

        # The last process gets any remaining particles.
        if process_index == number_processes - 1:
            end_index = len( unique_particle_ids )
        else:
            end_index = (process_index + 1) * particles_per_process

        batch_arguments.append( (unique_particle_ids[start_index:end_index],
                                 source_tree,
                                 destination_tree,
                                 keep_evaluations_flag) )

    # Keep track of the paths if we encounter I/O problems.
    read_error_paths  = []
    write_error_paths = []

    # Distribute the merge across the workers.
    with multiprocessing.Pool( processes=number_processes ) as pool:
        results = pool.starmap( merge_particle_batch, batch_arguments )

        for read_errors, write_errors in results:
            read_error_paths.extend( read_errors )
            write_error_paths.extend( write_errors )

    # Update destination particle index
    dest_index_path = droplet_approximation.get_particles_index_path( destination_tree )
    droplet_approximation.write_particles_index( dest_index_path,
                                                 unique_particle_ids,
                                                 merge_flag=True )

    # Report any errors that occurred
    if len( read_error_paths ) > 0:
        print( "Read errors occurred for the following files:\n" )
        for particle_path, exception_str in sorted( read_error_paths, key=lambda x: x[0] ):
            print( "  {:s} - {:s}".format( particle_path, exception_str ) )
        print()

    if len( write_error_paths ) > 0:
        print( "Write errors occurred for the following files:\n" )
        for particle_path, exception_str in sorted( write_error_paths, key=lambda x: x[0] ):
            print( "  {:s} - {:s}".format( particle_path, exception_str ) )
        print()

    if len( read_error_paths ) == 0 and len( write_error_paths ) == 0:
        print( "Directory tree merge completed successfully." )

def main( argv ):
    """
    Main entry point for the script.  Parses the command line and launches
    merge_directory_trees().

    Takes 1 argument:

      argv - Command line components, including the script's name as the first
             component.

    Returns 1 value:

      exit_code - Integral code to exit the script with.

    """

    usage_string = \
        "Usage: {:s} [-k] <source_tree> <destination_tree> <number_processes>".format(
            argv[0] )

    # Default to a safe position where evaluations are removed in the
    # destination directory.  This avoids the merged observations from being out
    # of sync and forces users to act unsafely rather than to act safely.
    keep_evaluations_flag = False

    # Parse any options (none currently used, but following the pattern)
    try:
        options, arguments = getopt.getopt( argv[1:], "k" )
    except getopt.GetoptError as error:
        print( "{:s}\n"
               "\n"
               "{:s}".format(
                   error,
                   usage_string ) )
        return 1

    for option, option_value in options:
        if option == "-k":
            keep_evaluations_flag = True
        else:
            raise NotImplementedError( "Unhandled option '{:s}'!".format( option ) )

    NUMBER_REQUIRED_ARGUMENTS = 3

    number_arguments = len( arguments )

    if number_arguments != NUMBER_REQUIRED_ARGUMENTS:
        print( usage_string )
        return 1

    source_tree      = arguments[0]
    destination_tree = arguments[1]
    number_processes = int( arguments[2] )

    # Validate the provided source tree path.
    if not os.path.exists( source_tree ):
        print( "Source tree ('{:s}') does not exist!".format( source_tree ) )
        return 1
    elif not os.path.isdir( source_tree ):
        print( "Source tree ('{:s}') is not a directory!}".format( source_tree ) )
        return 2

    # Create the top-level destination directory if it doesn't exist.
    os.makedirs( destination_tree, exist_ok=True )

    # Use one worker per core on the system if a specific worker count wasn't
    # provided.
    if number_processes == 0:
        number_processes = multiprocessing.cpu_count()

    # Do the merge.
    merge_directory_trees( source_tree,
                           destination_tree,
                           number_processes,
                           keep_evaluations_flag )

    return 0

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
