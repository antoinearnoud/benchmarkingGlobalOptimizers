#!/bin/bash
#SBATCH --job-name=localalgos_testfcts
#SBATCH --partition=day
#SBATCH --time=10:00:00
#SBATCH --output=job_output.txt
#SBATCH --cpus-per-task=25

# run optimization
./main

# copy results into results folder
cp results_local_algo_all.dat ../results/results_local_algo_all.dat
