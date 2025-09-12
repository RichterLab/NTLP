#!/usr/bin/env python3

# Script that reads a range of particles out of the supplied particles tree,
# possibly including evaluations, so as to verify consistency between the
# particles index and the trees contents.  Timing information is printed to
# standard output for diagnostics purposes (elapsed time, number of observations
# read, particles read rate, and overall read throughput).
#
# Usage: verify_particle_tree_consistency.py [-i <particles_index>] <particle_root> <particles_index_path> <start_particle_index> [<evaluation_tag> [<evaluation_tag> [...]]]

from datetime import datetime
import getopt
import sys
import time
from typing import List, Optional

import numpy as np

import droplet_approximation

# Number of directories per level in the particles hierarchy.  A branching
# factor of 256 provides roughly 65K directories to store particles in so as to
# not store too many particles in any one.  This is tuned for hierarchies
# containing millions of particles.
DIRS_PER_LEVEL = 256

def load_particles_dataframe( particles_root: str, particles_index_path: Optional[str],
                              particle_start_index: int, particle_end_index: int,
                              evaluation_tags_list: List[str] ):
    """
    Loads a subset of particles from the provided particles tree, possibly with
    their evaluation data, and reports timing information on the measured
    performance.  This is a simple function that serves as a canary in the coal
    mine to detect inconsistencies in a particles tree and the supplied index.
    Should an inconsistency be found, it is assumed this function will raise
    an exception identifying the fatal problem.

    Takes 5 arguments:

      particles_root       -
      particles_index_path -
      particle_start_index -
      particle_end_index   -
      evaluation_tags_list -

    Returns nothing.

    """

    # Get the path to the particles index we're using.
    if particles_index_path is None or len( particles_index_path ) == 0:
        particles_index_path = droplet_approximation.get_particles_index_path( particles_root )
    else:
        particles_index_path = "{:s}/{:s}".format(
            particles_root,
            particles_index_path )

    # Map the evaluation tags list into a evaluation tag map with dummy keys.
    # We don't actually care what the DataFrame column names are, only that
    # the evaluations associated with each tag are loaded correctly.
    evaluations = {}
    for evaluation_tag in evaluation_tags_list:
        evaluations[evaluation_tag] = evaluation_tag

    #
    # NOTE: The end index is exclusive so we don't need to add one.
    #
    number_particles = particle_end_index - particle_start_index

    # Report what we're doing and when.
    read_start_datetime = datetime.now()
    print( "Reading {:d} particle{:s} ({:d}:{:d}) from '{:s}' at {:s} (evaluations {:s}).".format(
        number_particles,
        "" if number_particles == 1 else "s",
        particle_start_index,
        particle_end_index,
        particles_index_path,
        read_start_datetime.strftime( "%Y/%m/%d %H:%M:%S" ),
        ", ".join( evaluation_tags_list ) ) )

    read_start = time.time()

    # Read the index and the particles requested to create a DataFrame.
    unique_particle_ids = np.fromfile( particles_index_path, np.int32 )
    particles_df        = droplet_approximation.read_particles_data(
        particles_root,
        unique_particle_ids[particle_start_index:particle_end_index],
        DIRS_PER_LEVEL,
        evaluations=evaluations )

    read_end          = time.time()
    read_end_datetime = datetime.now()

    total_number_evaluations = particles_df["number evaluations"].sum()
    total_number_seconds     = read_end - read_start

    # Report timing information
    print( "Finished at {:s} ({:.2f} seconds).  "
           "Read {:d} evaluation{:s} ({:.2f} particles/second, {:.1f} bytes/second).".format(
               read_end_datetime.strftime( "%Y/%m/%d %H:%M:%S" ),
               read_end - read_start,
               total_number_evaluations,
               "" if total_number_evaluations == 1 else "s",
               number_particles / total_number_seconds,
               total_number_evaluations * droplet_approximation.ParticleRecord.SIZE_BYTES / total_number_seconds ) )

def main( argv ):
    """
    Main entry point for the script.  Parses the command line and launches
    load_particles_dataframe().

    Takes 1 argument:

      argv - Command line components, including the script's name as the first
             component.

    Returns 1 value:

      exit_code - Integral code to exit the script with.


    """

    usage_string = \
        "Usage: {:s} [-i <particles_index>] <particles_root> <particles_start_index> <particles_end_index> [<evaluation_tag> [...]]".format(
            argv[0] )

    # Use the default particles index unless the user specifies on.
    particles_index_path = None

    # Only load the raw particle observations unless the user specifies
    # additional evaluations.
    evaluation_tags = []

    # Parse any options.
    try:
        options, arguments = getopt.getopt( argv[1:], "i:" )
    except getopt.GetoptError as error:
        print( "{:s}\n"
               "\n"
               "{:s}".format(
                   error,
                   usage_string ) )
        return 1

    for option, option_value in options:
        if option == "-i":
            particles_index_path = option_value
        else:
            raise NotImplementedError( "Unhandled option '{:s}'!".format( option ) )

    MINIMUM_NUMBER_ARGUMENTS = 3

    number_arguments = len( arguments )

    if number_arguments < MINIMUM_NUMBER_ARGUMENTS:
        print( usage_string )
        return 1

    # Parse the remaining arguments.
    particles_root       = arguments[0]
    particle_start_index = int( arguments[1] )
    particle_end_index   = int( arguments[2] )
    if number_arguments > MINIMUM_NUMBER_ARGUMENTS:
        evaluation_tags  = arguments[3:]

    # Check for consistency.
    load_particles_dataframe( particles_root,
                              particles_index_path,
                              particle_start_index,
                              particle_end_index,
                              evaluation_tags )

    return 0

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
