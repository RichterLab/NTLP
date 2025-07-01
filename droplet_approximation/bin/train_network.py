import functools
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

def generate_fortran_module_with_suffix( droplet_model_save_prefix, model_name, model, epoch_number ):
    """
    Serializes a model checkpoint's weights into a Fortran module suitable for
    inferencing.  This is a wrapper around generate_fortran_module() that is
    suitable to be called back during training so checkpoints have a companion
    Fortran module.

    Takes 4 arguments:

      droplet_model_save_prefix - Path prefix of the Fortran module to write.  This
                                  will have the epoch number appended to it, along
                                  with a ".f90" extension to construct the output
                                  path name.
      model_name                - Name of the model being serialized.
      model                     - PyTorch model being serialized.
      epoch_number              - Training epoch number associated with model.

    Returns nothing.

    """

    droplet_model_save_path = "{:s}_{:d}.f90".format(
        droplet_model_save_prefix,
        epoch_number )

    generate_fortran_module( droplet_model_save_path,
                             model_name,
                             model.state_dict() )

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
        epoch_callback = functools.partial( generate_fortran_module_with_suffix,
                                            droplet_model_save_prefix,
                                            model_name )
    else:
        epoch_callback = None

    # Train the model, potentially checkpointing and generating Fortran modules
    # along the way.
    training_loss = train_model( model,
                                 criterion,
                                 optimizer,
                                 device,
                                 number_epochs,
                                 training_data_path,
                                 model_save_prefix if epoch_checkpoint_flag else None,
                                 epoch_callback )

    print( "Trained on {:d} mini-batches.".format(
        len( training_loss ) ) )

    # Checkpoint the final model and generate a Fortran module for it.
    save_model_checkpoint( model_save_prefix,
                           -1,
                           model,
                           optimizer,
                           criterion,
                           training_loss )

    generate_fortran_module( droplet_model_save_path, model_name, model.state_dict() )
    print( "Wrote model weights and inferencing code to '{:s}'.".format( droplet_model_save_path ) )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
