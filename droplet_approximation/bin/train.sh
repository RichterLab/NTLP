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
#cleaning_script="./clean_data.py"
#tmp_data_dir="/scratch365/dcolange/pi_chamber/particle_traj/"
tmp_data_dir="../data/tmp"
data_file="../data/data_file_box_coupled_400M_stripped.data"
#clean_data_file="../data/data_file_NTLP_uncoupled_400M_stripped.data"
network_file="../models/coupled_epoch_checkpoints/network_box_coupled_400M_l1_full_residual.pth"
droplet_file="../models/coupled_epoch_checkpoints/droplet_model_box_coupled_400M_l1_full_residual.f90"

cat $tmp_data_dir/* > $data_file

#python3 ${cleaning_script} ${data_file} ${clean_data_file}
#python3 ${training_script} ${clean_data_file} ${network_file} ${droplet_file} -e
python3 ${training_script} ${data_file} ${network_file} ${droplet_file} -e

# Todo automate running data
