import warnings

import numpy as np
import torch
import torch.nn as nn
import os

import sys

from droplet_approximation.models import *

number_epochs = 14

def get_checkpoint_extension( model_save_path ):
    """
    Gets the file extension used in a checkpoint file's path.

    Takes 1 argument:

      model_save_path - Path to the model's checkpoint.

    Returns 1 value:

      extension_str - Extension on the checkpoint file's path.

    """

    return model_save_path.split( "." )[-1]

def main( argv ):
    if len( argv ) == 4:
        epoch_checkpoint_flag = False
    elif len( argv ) == 5 and argv[4] == "-e":
        epoch_checkpoint_flag = True
    else:
        print("Usage: <input_path> <network_output_path> <droplet_model_output_path> [<-e epoch checkpoint flag>]")
        sys.exit(1)

    training_data_path = argv[1]
    model_save_path = argv[2]
    droplet_model_save_path = argv[3]

    model_save_prefix = "{:s}_epoch".format( model_save_path.split( "." + get_checkpoint_extension( model_save_path ) )[0] )
    droplet_model_save_prefix = "{:s}_epoch".format( droplet_model_save_path.split( ".f90" )[0] )

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
    model     = ResidualNet()
    criterion = nn.L1Loss()
    #criterion = nn.MSELoss()
    optimizer = torch.optim.Adam( model.parameters(), lr=1.0e-3, weight_decay=1e-6 )
    
    # Move the model to the device we're training with.
    model = model.to( device )
    if epoch_checkpoint_flag:
        for i in range( number_epochs ):
            training_loss = train_model( model, 
                                        criterion,
                                        optimizer,
                                        device,
                                        1, 
                                        training_data_path )

            # Adjust learning rate
            for parameter_group in optimizer.param_groups:
                parameter_group["lr"] /= 2

            print( "Epoch {:n} trained on {:d} mini-batches.".format(
                i,
                len( training_loss ) ) ) 

            epoch_droplet_model_save_path = "{:s}_{:d}.f90".format(
                droplet_model_save_prefix,
                i + 1 )

            save_model_checkpoint( model_save_prefix,
                                   i,
                                   model,
                                   optimizer,
                                   criterion,
                                   training_loss )

            generate_fortran_module( model_name, model.state_dict(), epoch_droplet_model_save_path )
            print( "Wrote model weights and inferencing code to '{:s}'.".format( 
                epoch_droplet_model_save_path ) )

    else:
        training_loss = train_model( model, 
                                     criterion,
                                     optimizer,
                                     device,
                                     number_epochs, 
                                     training_data_path )
    
        print( "Trained on {:d} mini-batches.".format(
            len( training_loss ) ) ) 

        save_model_checkpoint( model_save_prefix,
                               -1,
                               model,
                               optimizer,
                               criterion,
                               training_loss )

        generate_fortran_module( model_name, model.state_dict(), droplet_model_save_path )
        print( "Wrote model weights and inferencing code to '{:s}'.".format( droplet_model_save_path ) )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
