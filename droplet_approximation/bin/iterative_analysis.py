
import torch
import functools

import sys

from droplet_approximation.models import *
from droplet_approximation.data import *
from droplet_approximation.analysis import parallel_analyze_model_iterative_performance_NTLP_data


def main( argv ):
    if len( argv ) != 5:
        print("Usage: <data_path> <model_load_path> <output_path> <iterations>")
        sys.exit(1)

    NTLP_data_path = argv[1]
    model_load_path = argv[2]
    output_path = argv[3]
    iterations = int(argv[4])

    device = torch.device( "cuda" if torch.cuda.is_available() else "cpu" )
    
    print( "Working with '{}'.".format( device ), flush=True)
    print("Loading model")
    
    model = SimpleNet()

    load_model_checkpoint( model_load_path, model )

    print("Loading data")
    data = pd.read_parquet(NTLP_data_path)

    print("Normalizing data")
    normalize_NTLP_df(data)
    
    print("Analyzing data")
    scaled_results = parallel_analyze_model_iterative_performance_NTLP_data( model, data, iterations, device, 64)

    print("Saving data")
    scaled_results.tofile(output_path)
    

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
