#!/usr/bin/env python3

# TODO:
#
#  1. Support MLP evaluation - determining the model class and providing a
#     model object to load_model_checkpoint() is required.
#  2. Correctly handle the MLP evaluation device.  "cuda" requires
#     multiprocessing to use spawn, though naively applying that here generates
#     a fork bomb.
#  3. Add a test mode to do everything but launch to let users verify before
#     creating a bunch of data they might not want.
#  4. Add a timer around the execution to report basic execution metrics.

import multiprocessing
import sys
import warnings

import numpy as np

import droplet_approximation

# Default to reviewing what would be done before we launch potentially
# hours/days worth of work.  Set this to False to execute as well.
testing_flag = True

# Manual toggle between BDF and MLP until this script is updated with proper
# command line flags.
mlp_flag = True

# Path to top-level directory in the NTLP repository.
repo_root = None

# Number of directories per level in the particle traces.
dirs_per_level = 256

if repo_root is None:
    raise RuntimeError( "Set data_root to a directory and run again!" )

if not mlp_flag:
    evaluation_type      = droplet_approximation.EvaluationType.BDF
    evaluation_extension = "bdf_iterative"
    iterative_flag       = True
    parameters           = {
        "atol":    (1e-10, 1e-4),
        "rtol":    (1e-7)
    }
else:
    evaluation_type      = droplet_approximation.EvaluationType.MLP
    evaluation_extension = "mlp-test"
    iterative_flag       = False
    parameters           = {
        "device":            "cpu",
        #
        # NOTE: This model wasn't saved as a v2 checkpoint so we have to provide
        #       its parameter ranges.
        #
        "checkpoint_path":   "{:s}/droplet_approximation/models/network_box_residual_l1_epoch_14.pt".format( repo_root ),
        "parameter_ranges":  {
            "radius": (-6.75, -3.00),
            "temperature": (284, 300),
            "salt_mass": (-17.66, -17.65),
            "air_temperature": (284.0, 300.0),
            "relative_humidity": (0.98, 1.11),
            "rhoa": (0.99, 1.01),
            "time": (-1.4, -0.75)
        }
    }

# Top-level path to the particles.
particles_root             = sys.argv[1]

# Indices within all available particles to process.
particle_indices_range_str = sys.argv[2]

# Number of processes to distribute the requested indices across.
number_processes           = int( sys.argv[3] )

particles_index_path = "{:s}/particles.index".format( particles_root )

# Default to all of the cores on the system if requested.
if number_processes == 0:
    number_processes = multiprocessing.cpu_count()

# Get the model's parameter ranges if we're doing MLP evaluation.
if evaluation_type == droplet_approximation.EvaluationType.MLP:
    # XXX: handle loading the model correctly
    model = droplet_approximation.ResidualNet()
    (parameter_ranges,
     _,
     _) = droplet_approximation.load_model_checkpoint( parameters["checkpoint_path"], model )

    parameters["model"]            = model

    if len( parameter_ranges ) == 0:
        warnings.warn( "'{:s}' did not contain parameter ranges - using defaults!".format(
            parameters["checkpoint_path"] ) )

        parameter_ranges = droplet_approximation.get_parameter_ranges()

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
particles_chunk_size = (len( particle_indices ) + number_processes - 1) // number_processes

args_list = []
particle_indices_list = []
for process_index in range( number_processes ):
    start_index = particles_chunk_size * process_index

    # Find the end of this process' range.  Take care to not exceed the last
    # requested particle while also handling any remaining particles if the
    # number of processes doesn't evenly divide the particle count.
    if process_index == number_processes - 1:
        end_index = min( len( particle_ids ),
                         particle_indices_range_components[1] )
    else:
        end_index = start_index + particles_chunk_size

    particle_indices_list.append( (start_index, end_index) )

    args = (particles_root,
            particle_ids[start_index:end_index],
            dirs_per_level,
            evaluation_type,
            evaluation_extension,
            iterative_flag,
            parameters)
    args_list.append( args )

# Report what we're processing.  Report particle numbers here.
print( "Processing particles {:d}-{:d} from '{:s}'.".format(
    particle_indices[0] + 1,
    particle_indices[-1],
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
else:
    raise ValueError( "Unknown evaluation type ({})!".format(
        evaluation_type ) )

# Stop if we're only interested in seeing what would be done.
if testing_flag:
    sys.exit( 0 )

# Launch the workers.
with multiprocessing.Pool( number_processes ) as pool:
    pool.starmap( droplet_approximation.particle_evaluation_pipeline, args_list )
