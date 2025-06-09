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

training_script="./continue_network_training.py"
data_file="../data/data_file_NTLP_decoupled_400M_stripped.data"
network_load_file="../models/network_NTLP_decoupled_400M_weighted.pth"
network_save_file="../models/network_NTLP_decoupled_400M_weighted_2.pth"
droplet_file="../models/droplet_model_NTLP_decoupled_400M_weighted_2.f90"

#python3 ${training_script} ${data_file} ${network_load_file} ${network_save_file} ${droplet_file} $SGE_TASK_ID
python3 ${training_script} ${data_file} ${network_load_file} ${network_save_file} ${droplet_file} 1

# Todo automate cleaning of data
