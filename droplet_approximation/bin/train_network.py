import warnings

import numpy as np
import torch
import torch.nn as nn
import os

import sys

from droplet_approximation.models import *

number_epochs = 14

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
            loss_history = train_model( model, 
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
                len( loss_history ) ) ) 

            epoch_model_save_path = "_epoch_{:n}.pth".format( i + 1 ).join(
                model_save_path.split(".pth") )
            epoch_droplet_model_save_path = "_epoch_{:n}.f90".format( i + 1 ).join(
                droplet_model_save_path.split(".f90") )

            torch.save( model.state_dict(), epoch_model_save_path )

            generate_fortran_module( model_name, model.state_dict(), epoch_droplet_model_save_path )
            print( "Wrote model weights and inferencing code to '{:s}'.".format( 
                epoch_droplet_model_save_path ) )

    else:
        loss_history = train_model( model, 
                                    criterion,
                                    optimizer,
                                    device,
                                    number_epochs, 
                                    training_data_path )
    
        print( "Trained on {:d} mini-batches.".format(
            len( loss_history ) ) ) 

        torch.save( model.state_dict(), model_save_path )

        generate_fortran_module( model_name, model.state_dict(), droplet_model_save_path )
        print( "Wrote model weights and inferencing code to '{:s}'.".format( droplet_model_save_path ) )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
