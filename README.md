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

All code for data generation, model training, and serialization can be found in in `ML/model_training`, including a checkpoint for our paper's model `models/flippant-gusto.pt`. `Model Training.ipynb` contains a full pipeline from data generation to model training that can be used to reproduce the model from our paper. The training/validation data used for this paper is present in `ML/model_training/data`. To reproduce the paper's model, be sure to train for the full `17,800` epochs rather than the default `100` used for testing.

Code for producing many of the figures used in the paper can be found at `ML/data_analysis/ten_history_post_process_comparison.ipynb`. This notebook requires that one download the paper data folder from the DOI and store it at `ML/data_analysis/paper_data`. 

Running the model on NTLP:
1. Copy the `droplet_model.f90` corresponding to your preferred network into the root directory.
2. Follow steps for compiling NTLP (quick note: run `make ARCH=avx2` when making): [link here](https://richterlab.miraheze.org/wiki/Setup_and_Running_NTLP)
3. Navitage to `test_cases/pi_chamber` and set `ipart_method=3` for neural network
4. If you want to dump the data from the standard NTLP simulation to generate training/testing data, keep `ipart_method=2` and set `iwritebe=1`. Control the fraction of data written with `iwritebe_proportion` as the traces can get very large very quickly.
5. Run `qsub pi_chamber.run`

For comparing BE/MLP results, you can use `ten_history_post_process_comparison.ipynb` by reducing `history_file_count` and pointing the notebook at the relavent BE/MLP history files.
