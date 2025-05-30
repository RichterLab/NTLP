# NTLP (NCAR Turbulence with Lagrangian Particles)

## Compilation
To build on the CRC machines you must run the following commands:
```
module load mvapich2
module load intel
module load netcdf

make clean
make
```

## SETUP AND RUNNING
To run, make a directory ("case1" or something) where les.run and params.in will go (i.e., not out of the same directory as the code files)

Set up the I/O directory on scratch: /scratch365/netID/tutorial/case1

Make sure all paths in params.in and les.run point to the proper locations


## For Droplet Approximation

The standard data generation and model training scripts are in droplet_approximation/bin. They are:
 * `generate_training_data.py` - Generates model training data. constants at the top `NUMBER_DROPLETS` and `NUMBER_SAMPLES` control this training process.
 * `parallel_generate_training_data.sh` - Runs an instance of generate_training_data.py for every core on the machine.
 * `train_network.py` - Usage: `python3 train_network.py <training_data_path> <model_save_path> <droplet_model_save_path>`. Trains a model for `number_epochs` epochs wiht data at `training_data_path`. Weights are recorded at `model_save_path` and the corresponding Fortran code is saved at `droplet_model_save_path`.
 * `train.sh` - Concatonates data from training data generation and runs `train_network.py` on said data.
 * `batch_queue.sh` - Usage: `./batch_queue <job_count>`. Queues `job_count` copies of `parallel_generate_training_data.sh` to the CRC cluster. Also queues an instance of `train.sh` to run after data generation is complete. Essentially automates a runthrough of data generation/model training on the CRC.

To generate data and train a model on the CRC:
1. Set desired parameter ranges in `droplet_approximation/python/droplet_approximation/physics.py`.
2. Set desired `NUMBER_SAMPLES` and `NUMBER_DROPLETS` in `approximation/generate_training_data.py`. Note that the overall amount of data generated will be `job_count X core_count X number_samples X number_droplets`. If you're not on the CRC presumably `job_count=1.`
3. Set desired output paths in `parallel_generate_training_data.sh` and `train.sh`.
4. If you're on the CRC, run `./batch_queue <job_count>`. It's worth playing a bit with how many jobs the queue will take at once. Often, it will do as many as 12. This makes data generation much faster.

Otherwise, train with `droplet_approximation/notebooks/Approximation Droplet Parameters.ipynb`. If desired, one can generate data on the CRC and export it for use with the notebook. 

Running on NTLP:
1. Copy the `droplet_model.f90` corresponding to your preferred network into the root directory.
2. Follow steps for compiling NTLP (quick note: run `make ARCH=avx2` when making): [link here](https://richterlab.miraheze.org/wiki/Setup_and_Running_NTLP)
3. Navitage to `test_cases/pi_chamber` and set `ipart_method=3` for neural network
4. If you want to dump the data from the standard NTLP simulation to generate training/testing data, keep `ipart_method=2` and set `iwritebe=1`. (BETA! Not yet working properly)
5. Run `qsub pi_chamber.run`

For visualizing results, use the matlab scripts provided with the results, namely `history_postprocessing.` To compare MLP output with standard backwards Euler, a reference history file `postprocessing/reference_history.nc` is provided with a script to compare the reference .nc file with the newly generated results: `postprocessing/compare_history_postprocessing`.

