#!/usr/bin/env python3

# Usage: generate_training_data.py <output_path> <spreadsheet_path>
#
# Creates 1,048,576 (1024^2) droplet parameters and writes them to <output_path>.
# Any weird parameters are logged to <spreadsheet_path> and ignored so they are
# not written to <output_path>.

import sys

from droplet_approximation import create_training_file

# XXX
NUMBER_DROPLETS = 4*1024 * 1024

# If 1, generate data normally.
# If > 1, will sample NUMBER_SAMPLES points logarithmically from a
#   8+DROPLET_TIME_LOG_RANGE(1) second integration of a droplet
#   Meant generation fo data sets with small time ranges.

NUMBER_SAMPLES = 1

def main( argv ):
    """
    """

    weird_file_name = None

    if len( argv ) == 1:
        print( "Usage: <output_file> [<weird_file>]" )
        sys.exit( 1 )
    if len( argv ) > 2:
        weird_file_name = argv[2]

    create_training_file( argv[1],
                          NUMBER_DROPLETS,
                          NUMBER_SAMPLES,
                          weird_file_name=weird_file_name )


if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
