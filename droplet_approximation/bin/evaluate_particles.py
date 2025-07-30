#!/usr/bin/env python3

# Usage: evaluate_particles.py [-f] <particles_root> <indices_range> <number_processes> <iterative_flag> bdf <extension> <atol_radii> <atol_temps> <rtol>
#                                   <particles_root> <indices_range> <number_processes> <iterative_flag> mlp <extension> <model_class> <device> <checkpoint_path>
#
# Evaluates a range of particles using N processes.  Evaluation is either BDF or
# MLP and each evaluation type's parameters are (painfully) specified on the
# commandline.  Since this has the potential for launching a lot of long-running
# work, it's default behavior is to report what would be done.  Execution is
# actually performed when the "-f" flag is provided.

# TODO:
#
#  1. Correctly handle the MLP evaluation device.  "cuda" requires
#     multiprocessing to use spawn, though naively applying that here generates
#     a fork bomb.
#  2. Add a timer around the execution to report basic execution metrics.

import getopt
import multiprocessing
import sys
import warnings

import numpy as np

import droplet_approximation

# Provide convenience variables with (slightly) shorter names.  We display these
# in lower-case for consistency with other arguments, regardless of their
# capitalization.
EVALUATION_TYPE_BDF_STR = droplet_approximation.EvaluationType.BDF.name.lower()
EVALUATION_TYPE_MLP_STR = droplet_approximation.EvaluationType.MLP.name.lower()

def create_new_model( model_class_name ):
    """
    Instantiates a PyTorch model object from a droplet_approximation's class
    name.  The object is created with default arguments.

    Raises ValueError if the supplied name does not exist in the
    droplet_approximation package or if it doesn't represent a class name.

    Takes 1 argument:

      model_class_name - String specifying the model class name to instantiate.
                         Must be one of the classes exposed in the
                         droplet_approximation package's interface.

    Returns 1 value:

      model - PyTorch model object described by model_class_name.

    """

    # Blow up if this isn't known to the package.
    if model_class_name not in droplet_approximation.__dict__:
        raise ValueError( "Model class '{:s} does not exist in droplet_approximation'!".format(
            model_class_name ) )

    # Lookup the symbol and confirm that it is indeed a class.
    obj = getattr( droplet_approximation, model_class_name )
    if not callable( obj ):
        raise ValueError( "Model class '{:s}' is not callable!".format(
            model_class_name ) )
    elif not isinstance( obj, type ):
        raise ValueError( "Model class '{:s}' is not a class!".format(
            model_class_name ) )

    # Instantiate it with the default arguments.
    model = obj()

    return model

def launch_evaluation_pipeline( particles_root, particle_indices_range_str,
                                number_processes, evaluation_type, iterative_flag,
                                evaluation_extension, testing_flag, parameters ):
    """
    Launches the particle evaluation pipeline using N-many processes.  This acts
    as a wrapper around droplet_approximation.particle_evaluation_pipeline()
    that converts command line arguments into the objects/values necessary for
    execution.

    Takes 8 arguments:

      particles_root             - Path to the particles/ sub-directory to process.
      particle_indices_range_str - String of the form <start>:<end> representing
                                   0-based particle indices.  <end> is exclusive
                                   to match Python range semantics, and may be
                                   -1 to represent the last available particle
                                   index.
      number_processes           - Non-negative integer specifying the number of
                                   processes to use during evaluation.  Zero
                                   requests one process per core on the system.
      evaluation_type            - droplet_approximation.EvaluationType
                                   enumeration specifying the type of
                                   evaluation.
      iterative_flag             - Boolean specifying whether iterative
                                   evaluation should be performed (True) or
                                   whether direct evaluation should be performed
                                   (False).
      evaluation_extension       - File extension to write the evaluation
                                   results to.  Existing files are overwritten.
      testing_flag               - Boolean specifying whether the configuration
                                   should only be printed (True) or if it should
                                   be printed and then executed (False).
      parameters                 - Dictionary of parameters for the requested
                                   evaluation type.

                                   BDF key/value pairs:

                                     atol:  Tuple of absolute tolerances, one for
                                            radii and one for temperatures.
                                     rtol:  Relative tolerance for both radii
                                            and temperatures.

                                   MLP key/value pairs:

                                     checkpoint_path:  Path to the model checkpoint path.
                                     device:           String specifying the device to
                                                       inference with.  Must be recognized
                                                       by PyTorch.
                                     model_class:      Name of the droplet_approximation
                                                       model class contained in checkpoint_path.

    Returns nothing.

    """

    # Number of directories per level in the particle traces.
    DIRS_PER_LEVEL = 256

    particles_index_path = "{:s}/particles.index".format( particles_root )

    # Default to all of the cores on the system if requested.
    if number_processes == 0:
        number_processes = multiprocessing.cpu_count()

    # Load the model and acquire its parameter ranges if we're doing
    # MLP evaluation.
    if evaluation_type == droplet_approximation.EvaluationType.MLP:
        model = create_new_model( parameters["model_class"] )

        # Promote warnings to errors so we can prevent v1 checkpoints from being
        # used.
        #
        # NOTE: While we could accept parameter ranges to facilitate evaluation
        #       with v1 checkpoints this is *very* error prone and we don't want
        #       to evaluate potentially billions of observations with the wrong
        #       values and generate total garbage.
        #
        with warnings.catch_warnings():
            warnings.simplefilter( "error", UserWarning )

            try:
                (parameter_ranges,
                 _,
                 _) = droplet_approximation.load_model_checkpoint( parameters["checkpoint_path"],
                                                                   model )
            except Warning as e:
                raise ValueError( "Checkpoint '{:s}' does not contain parameter ranges!\n"
                                  "Evaluation is only available for models that have explicit parameter ranges!".format(
                                      parameters["checkpoint_path"] ) )

        parameters["model"]            = model
        parameters["parameter_ranges"] = parameter_ranges

    # Figure out which particles are available.
    particle_ids = np.fromfile( particles_index_path, dtype=np.int32 )

    # Convert the range of particle indices into actual indices.
    particle_indices_range_components = list( map( lambda x: int( x ),
                                                   particle_indices_range_str.split( ":" ) ) )
    if len( particle_indices_range_components ) != 2:
        raise ValueError( "Only ranges with <start>:<stop> are supported ({:s})!".format(
            particle_indices_range_str ) )

    # Translate -1 to the actual last valid index.  Constructing a range() object
    # with an end of -1 creates an empty range and not through the end like a
    # slice().
    if particle_indices_range_components[1] == -1:
        particle_indices_range_components[1] = len( particle_ids )

        particle_indices = np.array( range( *particle_indices_range_components ) )

    if len( particle_indices ) == 0:
        raise ValueError( "Empty range specified! ({:s})".format(
            particle_indices_range_str ) )

    # Make sure the requested indices don't exceed the particles available.
    #
    # NOTE: We take the absolute value of the indices to handle the case
    #       where indexing from the end of the particles is requested.
    #
    maximum_index = np.abs( particle_indices ).max()
    if maximum_index >= len( particle_ids ):
        raise ValueError( "Indices requested are outside of particles available ([0:{:d}] < {:d})!".format(
            len( particle_ids ),
            maximum_index ) )

    # Break the requested indices up into process-sized chunks.
    # The last process gets the remainder.
    particles_chunk_size   = len( particle_indices ) // number_processes
    number_extra_particles = len( particle_indices ) % number_processes

    args_list             = []
    particle_indices_list = []
    start_index           = particle_indices[0]
    for process_index in range( number_processes ):
        end_index = (start_index + particles_chunk_size +
                     int( process_index < number_extra_particles ))

        # Check to see if we've run out of work for the remaining processes.
        # This occurs when there are fewer particles than processes.
        if start_index == end_index:
            print( "Only {:d} particle{:s} to process.  Reducing the number of processes from {:d} to {:d})".format(
                len( particle_indices ),
                "" if len( particle_indices ) == 1 else "s",
                number_processes,
                process_index  ) )
            number_processes = process_index

            break

        particle_indices_list.append( (start_index, end_index) )

        args = (particles_root,
                particle_ids[start_index:end_index],
                DIRS_PER_LEVEL,
                evaluation_type,
                evaluation_extension,
                iterative_flag,
                parameters)
        args_list.append( args )

        # Thanks to Python's exclusive indexing, the beginning of the next
        # processes' range is the end of this range.
        start_index = end_index

    # Report what we're processing.  Report particle numbers here.
    print( "Processing particles {:d}-{:d} from '{:s}'.".format(
        particle_indices[0] + 1,
        particle_indices[-1] + 1,
        particles_index_path ) )

    print( "\n"
           "          Particle Indices\n"
           "Process   Start      End     Count\n" )
    for process_index in range( number_processes ):
        print( "  {:3d}.  {:7d} - {:7d} {:8}".format(
            process_index + 1,
            particle_indices_list[process_index][0],
            particle_indices_list[process_index][1],
            # NOTE: These are indices so we don't need to add one.
            particle_indices_list[process_index][1] - particle_indices_list[process_index][0] ) )
    print()

    # Report the type of evaluation being performed.
    if evaluation_type == droplet_approximation.EvaluationType.BDF:
        print( "{:s} BDF evaluation:\n"
               "\n"
               "  Evaluation extension:  '{:s}'\n"
               "  Absolute tolerance:    {:s}\n"
               "  Relative tolerance:    {:s}\n".format(
                   "Iterative" if iterative_flag else "Direct",
                   evaluation_extension,
                   "SciPy default!" if "atol" not in parameters else "({:.1g}, {:.1g})".format( *parameters["atol"] ),
                   "SciPy default!" if "rtol" not in parameters else "{:.1g}".format( parameters["rtol"] )
               ) )
    elif evaluation_type == droplet_approximation.EvaluationType.MLP:
        print( "{:s} MLP evaluation:\n"
               "\n"
               "  Evaluation extension:  '{:s}'\n"
               "  Model path:            {:s}\n"
               "  Inference device:      {:s}\n".format(
                   "Iterative" if iterative_flag else "Direct",
                   evaluation_extension,
                   parameters["checkpoint_path"],
                   parameters["device"],
               ) )
        droplet_approximation.display_parameter_ranges( parameters["parameter_ranges"],
                                                        droplet_approximation.DisplayType.HUMAN,
                                                        "  " )
        print()

    # Stop if we're only interested in seeing what would be done.
    if testing_flag:
        return 0

    # Otherwise, launch the workers.
    with multiprocessing.Pool( number_processes ) as pool:
        pool.starmap( droplet_approximation.particle_evaluation_pipeline,
                      args_list )

def main( argv ):
    """
    """

    testing_flag = True

    #
    # NOTE: We track BDF and MLP separately to handle the (likely) event that
    #       they will require different numbers of arguments.
    #
    NUMBER_BDF_ARGUMENTS = 9
    NUMBER_MLP_ARGUMENTS = 9

    number_arguments = len( argv[1:] )

    # We have two different calling interfaces, one for BDF and one for MLP.
    # We print them on separate lines instead of trying to show them
    # on one and confusing everyone.
    usage_string = \
        "Usage: {:s} [-f] <particles_root> <indices_range> <number_processes> <iterative_flag> bdf <extension> <atol_radii> <atol_temps> <rtol>\n" \
        "       {:>{:d}} [-f] <particles_root> <indices_range> <number_processes> <iterative_flag> mlp <extension> <model_class> <device> <checkpoint_path>".format(
            argv[0],
            "",
            len( argv[0] ) )

    if (number_arguments != NUMBER_BDF_ARGUMENTS) and (number_arguments != NUMBER_MLP_ARGUMENTS):
        print( usage_string )
        return 1

    # Parse any options to handle the case where we're "forced" to execute.
    try:
        options, arguments = getopt.getopt( sys.argv[1:], "f" )
    except getopt.GetoptError as error:
        print( "{:s}\n"
               "\n"
               "{:s}".format(
                   error,
                   usage_string ) )
        return 1

    for option, option_value in options:
        if opt == "-f":
            testing_flag = False
        else:
            raise NotImplementedError( "Unhandled option '{:s}'!".format( option ) )

    particles_root   = arguments[0]
    indices_range    = arguments[1]
    number_processes = int( arguments[2] )
    iterative_flag   = bool( int( arguments[3] ) )
    evaluation_type  = arguments[4]
    file_extension   = arguments[5]


    # Convert the command line evaluation type string into the internal
    # constants and package each type's parameters into a dictionary.
    parameters = {}
    if evaluation_type.lower() == EVALUATION_TYPE_BDF_STR:
        evaluation_type    = droplet_approximation.EvaluationType.BDF
        parameters["atol"] = (float( arguments[6] ), float( arguments[7] ))
        parameters["rtol"] = float( arguments[8] )
    elif evaluation_type.lower() == EVALUATION_TYPE_MLP_STR:
        evaluation_type               = droplet_approximation.EvaluationType.MLP
        parameters["model_class"]     = arguments[6]
        parameters["device"]          = arguments[7]
        parameters["checkpoint_path"] = arguments[8]
    else:
        raise ValueError( "Unknown evaluation type ('{:s}')!  Must be either '{:s}' or '{:s}'.".format(
            evaluation_type,
            droplet_approximation.EvaluationType.BDF.name.lower(),
            droplet_approximation.EvaluationType.MLP.name.lower() ) )

    # Execute the pipeline.
    launch_evaluation_pipeline( particles_root,
                                indices_range,
                                number_processes,
                                evaluation_type,
                                iterative_flag,
                                file_extension,
                                testing_flag,
                                parameters )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
