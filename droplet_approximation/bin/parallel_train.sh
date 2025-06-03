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

index=$SGE_TASK_ID

batch_name="NTLP-pull"

training_script="./train_network.py"
data_file="../data/data-file.data"
network_file="../models/network-${batch_name}-${index}.pth"
droplet_file="../models/droplet_model-${batch_name}-${index}.f90"

#rm $tmp_data_dir/*.dat
#mv $tmp_data_dir/-p $tmp_data_dir/../tmp_unsure
#cat $tmp_data_dir/* > $data_file

python3 ${training_script} ${data_file} ${network_file} ${droplet_file}
