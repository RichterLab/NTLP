#!/usr/bin/env python3

# Usage: generate_training_data.py <output_prefix> [<number_processes> [<number_droplets> [<seed> [<ranges>]]]]
#
# Creates random droplet data in parallel using <number_processes>-many
# processes.  Each process creates <number_droplets>-many droplets using the RNG
# seed <seed> (offset by the process' rank), with each droplet sampled from the
# <ranges> distributions.  Everything but the output prefix is optional and
# provides defaults.

import multiprocessing
import sys

import numpy as np

import droplet_approximation

# ~1M droplets is a sizable chunk of work to perform.
DEFAULT_NUMBER_DROPLETS = 1024 * 1024

# Default to serial execution.
DEFAULT_NUMBER_PROCESSES = 1

# Default to each process setting a random seed.
#
# NOTE: This results in streams of non-reproducible random values!
#
DEFAULT_RNG_SEED = None

# Use the droplet_approximation package's parameter ranges.
DEFAULT_PARAMETER_RANGES = None

# Names of each of the parameter ranges, in the order they must
# be provided on the command line.
PARAMETER_NAMES = ["radius",
                   "temperature",
                   "salt_solute",
                   "air_temperature",
                   "relative_humidity",
                   "rhoa",
                   "time"]

def parse_parameter_ranges( parameter_ranges_str ):
    """
    Parses a string containing the droplet_approximation package's parameter
    ranges and returns a dictionary containing the parsed bounds.  Each
    parameter's range is colon-delimited, with each range separated by commas
    like so:

      <radius_low>:<radius_high>,<temp_low>:<temp_high>,<solute_low>:<solute_high>,<air_temp_low>:<air_temp_high>,<rh_low>:<rh_high>,<rhoa_low>:<rhoa_high>,<time_low>:<time_high>

    All seven of the parameters must have ranges specified.

    Takes 1 argument:

      parameter_ranges_str - Parameter ranges string to parse.  See above for the
                             format.

    Returns 1 value:

      parameter_ranges - Dictionary of parameter names and a sequence of length
                         two containing their lower and upper bounds.  Parameter
                         names must be compatible with droplet_approximation.set_parameter_ranges().

    """

    if parameter_ranges_str is None:
        return {}

    parameter_components = parameter_ranges_str.split( "," )

    # All of the parameters must have ranges specified.  This forces the caller
    # to make a conscious decision about the ranges instead of assuming some
    # defaults and overriding others.
    if len( parameter_components ) != 7:
        raise ValueError( "Expected the parameter ranges to have 7 components, received {:d}.".format(
            len( parameter_components ) ) )

    parameter_ranges = {}
    for parameter_index, parameter_name in enumerate( PARAMETER_NAMES ):
        # Get the bounds as floating point values.
        parameter_range = list( map( lambda parameter: float( parameter ),
                                     parameter_components[parameter_index].split( ":" ) ) )

        # Ensure we only have a lower and upper bound.
        if len( parameter_range ) != 2:
            raise ValueError( "Expected '{:s}' to have 2 components (start, stop) but received {:d}.".format(
                parameter_components[parameter_index],
                len( parameter_range ) ) )

        # Track this range.
        parameter_ranges[parameter_name] = parameter_range

    return parameter_ranges

def create_training_file_with_seed_and_ranges( output_path, log_path, number_droplets, seed, parameter_ranges ):
    """
    Sets the NumPy RNG seed and the droplet parameter ranges and creates a file
    of random droplets.  After completion a summary message is printed to
    standard output.  Standard error is logged to file so failed ODE solves can
    be analyzed.

    Takes 5 arguments:

      output_path      - Path to the file to create.  This is overwritten if it
                         already exists.
      log_path         - Path to the log file to create.  Standard output is
                         redirected here.
      number_droplets  - The number of droplets to generate.
      seed             - Integral seed for NumPy's RNG.
      parameter_ranges - Dictionary of parameter ranges for to set with
                         set_parameter_ranges().

    Returns nothing.

    """

    # Set a global seed to what we were provided.
    #
    # NOTE: This isn't ideal, or the preferred approach, but is the best we've
    #       got to influence functions that generate random numbers without
    #       exposing the RNG as part of their interface.
    #
    np.random.seed( seed )

    droplet_approximation.set_parameter_ranges( parameter_ranges )

    # Redirect standard error to the log and create the training file.
    with open( log_path, "w" ) as log_fp:
        sys.stderr = log_fp

        droplet_approximation.create_training_file( output_path,
                                                    number_droplets,
                                                    quiet_flag=False )

    print( "Wrote '{:s}' with {:d} droplet{:s}.".format(
        output_path,
        number_droplets,
        "" if number_droplets == 1 else "s" ) )

def main( argv ):
    """
    Generates random droplet data for training, validation, and testing.

    Takes 1 argument:

      argv - Sequence of command line arguments to parse, whose first argument is
             the path to the script.  Must have between 2 and 6 elements:

               1. Script path name
               2. Output prefix for data generated
               3. Optional number of processes to generate data in parallel.
                  Must be a non-negative integer, zero specifies use one process
                  per core on the system.  If omitted, defaults to
                  DEFAULT_NUMBER_PROCESSES.
               4. Optional number of droplets to generate per process.  If
                  omitted, defaults to DEFAULT_NUMBER_DROPLETS.
               5. Optional RNG seed.  If omitted, defaults to DEFAULT_RNG_SEED.
               6. Optional parameter range string.  If omitted, defaults to
                  DEFAULT_PARAMETER_RANGES.

    Returns 1 value:

      exit_status - Integral exit status.  Zero is success, non-zero indicates
                    a failure of some kind.

    """

    NUMBER_MIN_ARGUMENTS = 2
    NUMBER_MAX_ARGUMENTS = 6

    number_arguments = len( argv )

    if (NUMBER_MIN_ARGUMENTS > number_arguments) or (NUMBER_MAX_ARGUMENTS < number_arguments):
        print( "Usage: {:s} <output_prefix> [<number_processes> [<number_droplets> [<seed> [<ranges>]]]]".format(
            argv[0] ) )

        return 1

    # Parse each of the arguments and default the ones not provided.
    output_prefix = argv[1]
    if number_arguments > 2:
        number_processes = int( argv[2] )
    else:
        number_processes = DEFAULT_NUMBER_PROCESSES
    if number_arguments > 3:
        number_droplets = int( argv[3] )
    else:
        number_droplets = DEFAULT_NUMBER_DROPLETS
    if number_arguments > 4:
        rng_seed = int( argv[4] )
    else:
        rng_seed = DEFAULT_RNG_SEED
    if number_arguments > 5:
        parameter_ranges_str = argv[5]
    else:
        parameter_ranges_str = DEFAULT_PARAMETER_RANGES

    # Use one process per core available if requested.
    if number_processes == 0:
        number_processes = multiprocessing.cpu_count()

    # Overlay the requested parameter ranges on top of the defaults.
    parameter_ranges = droplet_approximation.get_parameter_ranges()
    parameter_ranges.update( parse_parameter_ranges( parameter_ranges_str ) )

    print( "Generating {:d} file{:s} of data with the prefix '{:s}' each with {:d} droplet{:s} ({:d} total).\n".format(
        number_processes,
        "" if number_processes == 1 else "s",
        output_prefix,
        number_droplets,
        "" if number_droplets == 1 else "s",
        number_processes * number_droplets ) )

    # Print out the ranges so the user can review them.
    droplet_approximation.display_parameter_ranges( parameter_ranges )
    print()

    # Spawn off N-many processes each generating a block of data.
    with multiprocessing.Pool( processes=number_processes ) as pool:
        args_list = []

        for process_index in range( number_processes ):

            # Build this process' pathes.
            output_path = "{:s}-{:d}.training_data".format( output_prefix, process_index )
            log_path    = "{:s}-{:d}.log".format( output_prefix, process_index )

            # Adjust the provided seed by this process' rank so it has a
            # distinct stream of values from its brethren.  If no seed is
            # provided, each process will seed its NumPy RNG from random data.
            if rng_seed is not None:
                current_rng_seed = rng_seed + process_index
            else:
                current_rng_seed = None

            args_list.append( (output_path,
                               log_path,
                               number_droplets,
                               current_rng_seed,
                               parameter_ranges ) )

        pool.starmap( create_training_file_with_seed_and_ranges,
                      args_list )

    return 0


if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
