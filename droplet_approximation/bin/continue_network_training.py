#!/usr/bin/env python3

import warnings

import numpy as np
import torch
import torch.nn as nn
import os

import sys

from droplet_approximation.models import *

configs = [[1,1.0e-8, True]]
#        [3,0.2e-7, True],
#        [3,0.4e-8, False],
#        [3,0.4e-8, True]]
#

def main( argv ):
    if len( argv ) != 6:
        print("Usage: <data_input_path> <network_input_path> <network_output_path> <droplet_model_output_path> <config>")
        sys.exit(1)

    config_index = int(argv[5]) - 1
    training_data_path = argv[1]
    model_load_path = argv[2]
    model_save_path = argv[3] + str(config_index)
    droplet_model_save_path = argv[4] + str(config_index)

    config = configs[config_index]

    number_epochs = config[0]

    model_name = model_save_path.split("/")[-1].split(".")[0]

    device = torch.device( "cuda" if torch.cuda.is_available() else "cpu" )
    
    print( "Training with '{}'.".format( device ) )
    
    # Create an instance of SimpleNet and configure its optimization parameters:
    #
    #   - We use mean-squared error loss so we chase outliers
    #   - We use Adam so we have momentum when performing gradient descent
    #   - Given that we have a relatively large batch size we use a large
    #     initial learning rate of 1e-3.
    #   - We want smaller weights so we use a L2 regularization penalty of 1e-6
    #     to encourage that (as well as having non-zero weights rather than
    #     leaning heavily on just a subset of weights)
    #
    model     = SimpleNet()
    criterion = nn.MSELoss()
    #criterion = weighted_mse_loss
    if config[2]:
        optimizer = torch.optim.Adam( model.parameters(), lr=config[1], weight_decay=1.0e-6)
    else:
        optimizer = torch.optim.Adam( model.parameters(), lr=config[1])

    load_model_checkpoint( model_load_path, model, optimizer )
    
    # Move the model to the device we're training with.
    model = model.to( device )
    training_loss = train_model( model,
                                 criterion,
                                 optimizer,
                                 device,
                                 number_epochs,
                                 training_data_path )

    print( "Trained on {:d} mini-batches.".format(
        len( training_loss ) ) )

    torch.save( model_save_path,
                -1,
                model,
                optimizer,
                criterion,
                training_loss )

    generate_fortran_module( droplet_model_save_path, model_name, model.state_dict() )
    print("Wrote model weights and inferencing code to '{:s}'.".format(droplet_model_save_path))

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
