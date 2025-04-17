# Overview
This directory contains training data, though none is included as part of the
repository due to its size (O(10 GiB).  Data can be generated, directly copied
or linked to as necessary.  The model training process reads the entirety of the
data set before the first epoch so pick the highest I/O performance approach
possible given your storage limitations (if any).

# Creating Training Data
The Jupyter Notebook's `create_training_file()` can be used to create new
training data.  As of October of 2024 this method is single-threaded though
there is no technical reason why it could not be run in parallel using
the `multiprocessing` module's process pools or a distributed task system like
Dask.

In fact the `bin/parallel_generate_training_data.sh` provides this as a
standalone script and generates a number of files in batches, one per core on
the local system.  Generating data on a single core can also be accomplished
with `bin/generate_training_data.py`.  Until a better way of managing the
training workspace, Both of these should be run from inside the `bin/`
sub-directory.

``` shell
$ cd ../bin/
$ ./parallel_generate_training_data.sh
Generating log data #0000
Generating log data #0001
Generating log data #0002
Generating log data #0003
Generating log data #0004
Generating log data #0005
Generating log data #0006
Generating log data #0007
Generating log data #0008
Generating log data #0009
Generating log data #0010
Generating log data #0011
Generating log data #0012
Generating log data #0013
Generating log data #0014
Generating log data #0015
...
$ ls ../data/time_log_spaced-*.data
time_log_spaced-0000.data  time_log_spaced-0001.data  time_log_spaced-0003.data
time_log_spaced-0004.data  time_log_spaced-0005.data  time_log_spaced-0006.data
...
```

## Combining Training Data
The training process could be updated to take a set of files that are
individually read and concatenated in memory but this requires additional
time that was not available during initial development.

Data generated with any of these methods will create one or more files
containing droplet parameters.  These must be combined into a single file
before training which is luckily as easy as using the `cat` command:

``` shell
# NOTE: This overwrites time_log_spaced-combined.data!
$ cat ../data/time_log_spaced-*.data ... > ../data/time_log_spaced.data
```

## Weird Droplet Parameters
The ODEs governing droplets size and temperature may generate data that are
outside of the "reasonable" parameter space.  Some of these seem plausible
given their inputs though said inputs are physically unrealistic (either
outright impossible or outside the simulations of interest) while others
are physically impossible (e.g. negative radii or temperature).

Whenever these are encountered they are logged and written out as an Excel
spreadsheet so they can be analyzed and the underlying ODEs refined.  By default
said spreadsheets are generated next to the training data, one per file.
