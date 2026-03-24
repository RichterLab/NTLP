#!/bin/bash

#$ -M dcolange@nd.edu
#$ -m n
#$ -pe smp 64
#$ -q *@@richter
#$ -N optimize_tvm
#$ -j y
#$ -e ./logs/
#$ -o ./logs/

# make sure relavent libraries are installed (torch, matplotlib, scipy, numpy, etc)

conda activate NTLP

python3 optimize_tvm.py
