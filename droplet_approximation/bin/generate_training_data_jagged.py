#!/usr/bin/env python3

# Usage: generate_training_data.py <output_path> <spreadsheet_path>
#
# Creates 1,048,576 (1024^2) droplet parameters and writes them to <output_path>.
# Any weird parameters are logged to <spreadsheet_path> and ignored so they are
# not written to <output_path>.

import sys

from droplet_approximation import create_training_file

# XXX
NUMBER_DROPLETS = 512*1024
NUMBER_SAMPLES = 2

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
