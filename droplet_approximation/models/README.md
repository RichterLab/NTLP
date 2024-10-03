# Overview
This directory contains trained models that should be retained for one reason
or another.

Below are details on previously trained models, the design decisions that went
into their architecture and hyperparameter selection, the training data they
used, and how they were serialized for retention and future use.

***NOTE:*** Please be mindful of checking in new models as they take up space
in the repository that is not reclaimed should the model be removed after the
fact.  For the small models initially considered, this is not a serious concern
but this warning remains as a reminder should larger models be explored in the
future.

# Models
Of all of the models initially trained, two have been retained with the
following nicknames:

1. [Batch=256, LR=1e-3](#model-1)
2. [Batch=1024, LR=1e-3+schedule, L2=1e-6](#model-2)

Both of these have good performance at estimating the underlying ODEs they're
trained on.  No quantitative comparison has been performed but qualitatively
the Batch=256, LR=1e-3 model has better agreement of the size and temperature
outputs across the target time interval [1e-3, 10] seconds, though when it
estimates the wrong values it is "more" wrong than the other model.

See the corresponding sections for additional details on how they were trained.

## Design Choices
Very little energy has been spent selecting a model architecture other than
using one that will execute quickly at run-time.  If estimation of droplet
parameters (model inference) is slower than the previous Gauss-Newton iterative
approach then there is no point in using a neural network as it a) slows down
the overall simulation and b) produces less accurate droplet parameters.

The initial selection of a 4-layer multilayer perceptron (MLP) with ReLU
activations, with each layer having 32 neurons, was an attempt to have a small
network as well as overall learnable parameter count.  A small network impacts
inference efficiency as each layer corresponds to a matrix multiplication and
parameter count has implications on the accuracy of the solution.  Based on
initial (flawed) timings, this was faster than the Gauss-Newton iterative
approach on both Intel Ivy Bridge and AMD Epyc nodes and effort went towards
training the model to a sufficient level of accuracy.

As training data generation was refined to only include droplet parameters
that produced physically possible outputs (e.g. ignoring negative radii and
temperatures) a very small hyperparameter search was performed.  This was
motivated by the fact that the [best model produced at the
time](#model-1)
could not be reproduced with subsequent training runs.  As a result the
following areas were briefly explored:

- Increasing the batch size from 256
- Reducing the learning rate
- Adding a learning rate schedule
- Adding weight regularlization to the loss function

***NOTE:*** There was insufficient time to rigorously evaluate hyperparameter
selection performance and it was done in an ad hoc manner without quantiatively
measuring a model's performance.  As such, take the following as a statement of
what was done rather than one of fact.

Given that the network, its inputs and outputs are all small the batch size was
increased from 256 to see if less gradient variance would improve the quality.
Batches sized 512, 1024, 2048, 5120, and 10240 were explored with various
learning rates but, qualitatively, performance was poorer when batches were
larger than 1024.

A learning rate of 1e-3 is relatively large though appears to work reasonably
well.  Reducing the learning rate, with 5e-3, 1e-4, 1e-5, and 1e-6, produced
well behaved models but did not approach the accuracy of a 1e-3 learning rate
model.

While a learning rate of 1e-3 generally produced the best models it also
resulted in quite a bit of variance in the training loss.  As such, a halving
schedule was introduced so that the learning rate was cut in half for each
epoch of training.  For a 10 epoch training run this admitted a large learning
rate to get weights into the correct area with big initial steps and refined
them with smaller steps in subsequent epochs.

An attempt to improve model performance by encouraging small weights introduced
a L2 regularization penalty into the loss equation.  This was introduced around
the same time as the learning rate schedule and was not assessed independently.
It does not appear to have hurt performance.

## Model #1
This was the second model ever trained.  It is a 4-layer MLP with ReLU
activations with the following layer input sizes:

Input: 7
Layer 1: 32
Layer 2: 32
Layer 3: 32
Output: 2

This model has 2432 trainable parameters.  It resides at `models/mlp_4layer-b=256-lr=1e-3.pth`.

It was trained while droplet parameters were still being sampled and its
training set only had 18.3 million droplets in it.  The t_final parameter
was [sampled in log-space as described below](#Training_Data).  The model
was trained on 10 epochs of this dataset, with each epoch shuffling the
droplet parameters' order.  Each epoch was processed with batches of 256
droplet parameters.

A relatively large learning rate of 1e-3 was fixed throughout the training
run.

The model weights at the end of 10 epochs were serialized to disk.  No
checkpoints were retained during the training so there was no ability to select
a model with the best loss within a region at the end of training.

## Model #2
This was one of the last models trained.  It is a 4-layer MLP with ReLU
activations with the following layer input sizes:

Input: 7
Layer 1: 32
Layer 2: 32
Layer 3: 32
Output: 2

This model has 2432 trainable parameters.  It resides at `models/mlp_4layer-b=1024-lr=1e-3_halvingschedule-l2reg=1e-6.pth`.

It was trained after droplet parameters were sampled and had a training
set of roughly 309.9 million droplets in it.  The t_final parameter
was [sampled in log-space as described below](#Training_Data).  The model
was trained on 10 epochs of this dataset, with each epoch shuffling the
droplet parameters' order.  Each epoch was processed with batches of 1024
droplet parameters.

A relatively large initial learning rate of 1e-3 was used.  Each epoch had
the learning rate halved allowing earlier epochs to make larger weight updates
while later epochs were constrained to much smaller refinements (O(1000x)
smaller in the 10th epoch).

The model weights at the end of 10 epochs were serialized to disk.  No
checkpoints were retained during the training so there was no ability to select
a model with the best loss within a region at the end of training.

# Training Data
Training data were generated sampling the interval [-1, 1] and mapping it into
physical units for the non-temporal droplet parameters (radius, temperature,
salinity, air temperature, relative humidity, and rhoa).  The estimation time,
t_final, was sampled in log-space over the exponent interval of [-3.2, 1.1].
While the temporal ranges anticipated in DNS and LES span the exponent interval
of [-3, 1] the extra range was used to learn the dynamics at the interval
endpoints so it would not extrapolate (as much) when estimating the range
extremes.

Data were generated in parallel, on multiple cores to produce roughly 335.5
million random droplets that were concatenated together a single file.  Note
that roughly 25 million of said parameters were deemed unsuitable for training
(e.g. physically impossible or well outside of the range of physical parameters)
so only 309.9 million droplet parameters should be used.  See the
`clean_training_data()` function for details on how droplet parameters were
filtered.

# Formats
No specific format needed, though consider the dependencies required for working
with a model.

Serialized PyTorch models (.pth files) are .zip files containing the model
weights and metadata describing the layers they're associated with.  These can
be saved and loaded with just PyTorch and no additional dependencies.

Compare this to models serialized in the ONNX file format which requires the
ONNX software stack in addition to PyTorch.  This may be acceptable in the
future should the models need to be used in a variety of tools but was avoided
due to the limited scope of use during initial development (e.g. train a model
and serialize its weights into a Fortran module for inference, maybe analyze
said model's performance within a Jupyter Notebook).
