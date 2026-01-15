#!/usr/bin/env python3

# Usage: score_report.py config_path simulation_profile error_profile
#
# Generates a score report for the given simulation with the settings
# specified in the error config.
#
# The resulting score report will be stored in a folder stored at the location
# specified by the "output_root" option in the error profile.

import sys
import pickle
import os

# Limit numpy to 1 thread so that
# we can parallelize the error analysis
# properly

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import droplet_approximation

def main( argv ):
    if len( argv ) != 4:
        print( "Usage: <config_path> <simulation_profile> <error_profile>" )

    # Set up and load config
    config_path        = argv[1]
    simulation_profile = argv[2]
    error_profile      = argv[3]

    config = droplet_approximation.ExtendedConfigParser()
    config.load_config( config_path )
    config.load_subconfig( "simulation", simulation_profile )
    config.load_subconfig( "error_analysis", error_profile )

    error_config      = config["error_analysis"]
    simulation_config = config["simulation"]

    print( "Running Error Analysis with the following config:" )
    config.display()

    # Create output directories
    reference_tag  = error_config.get( "reference_tag" )
    comparison_tag = error_config.get( "comparison_tag" )

    analysis_root_dir = config["error_analysis_output_paths"].get( "output_root" ) + "/{:s}/".format( comparison_tag )
    analysis_dir_name = droplet_approximation.generate_analysis_directory_name( config )

    analysis_dir_path          = analysis_root_dir + analysis_dir_name
    figures_dir_path           = analysis_dir_path + "figures/"
    deviation_figures_dir_path = figures_dir_path  + "deviations/"
    score_report_path          = analysis_dir_path + "score_report.pkl"
    substitutions_dir_path     = analysis_dir_path + "substitutions/"

    os.makedirs(figures_dir_path, exist_ok=True)
    os.makedirs(substitutions_dir_path, exist_ok=True)
    os.makedirs(deviation_figures_dir_path, exist_ok=True)

    # Run score report and pickle the result
    score_report = droplet_approximation.ScoreReport( config )

    with open( score_report_path, "wb" ) as score_file:
        pickle.dump( score_report, score_file )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
