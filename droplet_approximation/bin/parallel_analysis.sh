#!/bin/bash

#$ -M dcolange@nd.edu
#$ -m n 
#$ -pe mpi-64 64
#$ -q *@@richter
#$ -N analysis 
#$ -j y
#$ -e ./logs/
#$ -o ./logs/

conda activate NTLP

data_path="../data/analysis_df.parquet"
output_path="../data/error_score_NTLP_data.data"
model_path="../models/network-NTLP.pth"
iterations = 5

analyzer="./iterative_analysis.py"

python3 -u ${analyzer} ${data_path} ${model_path} ${output_path} ${iterations}
