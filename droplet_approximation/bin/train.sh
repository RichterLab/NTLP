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

module unload python

training_script="./train_network.py"
tmp_data_dir="../data/tmp"
data_file="../data/data-file.data"
network_file="../models/network.pth"
droplet_file="../models/droplet_model.f90"

cat $tmp_data_dir/* > $data_file
rm $tmp_data_dir -R

python3 ${training_script} ${data_file} ${network_file} ${droplet_file}
