#!/bin/bash

# Usage: score_report.sh config_path simulation_profile error_profile
#
# Generates a score report for the given simulation with the settings
# specified in the error config. 
#
# The resulting score report will be stored in a folder stored at the location
# specified by the "output_root" option in the error profile.

#$ -M dcolange@nd.edu
#$ -m n 
#$ -pe mpi-64 64
#$ -q *@@richter
#$ -j y
#$ -e ./logs/
#$ -o ./logs/

conda activate NTLP

python3 score_report.sh $1 $2 $3
