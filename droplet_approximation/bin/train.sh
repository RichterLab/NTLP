#!/bin/bash

#$ -M dcolange@nd.edu
#$ -m n 
#$ -pe mpi-64 64
#$ -q *@@richter
#$ -N train_network 
#$ -j y
#$ -e ./logs/
#$ -o ./logs/

# make sure relavent libraries are installed (torch, matplotlib, scipy, numpy, etc)

conda activate NTLP

training_script="./train_network.py"
cleaning_script="./clean_data.py"
tmp_data_dir="/scratch365/dcolange/pi_chamber/particle_traj/"
data_file="../data/data_file_NTLP_uncoupled_400M.parquet"
clean_data_file="../data/data_file_NTLP_uncoupled_400M_stripped.data"
network_file="../models/uncoupled_mse_checkpoints/network_NTLP_uncoupled_400M.pth"
droplet_file="../models/uncoupled_mse_checkpoints/droplet_model_NTLP_uncoupled_400M.f90"

#rm $tmp_data_dir/*.dat
#cat $tmp_data_dir/* > $data_file

#python3 ${cleaning_script} ${data_file} ${clean_data_file}
python3 ${training_script} ${clean_data_file} ${network_file} ${droplet_file} -e

cp ${droplet_file} ../../droplet_model.f90
cd ../../

make ARCH=avx2

cd test_cases/pi_chamber

qsub pi_chamber.run

# Todo automate cleaning of data
