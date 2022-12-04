#!/bin/bash
#SBATCH --job-name=nlopt_testfcts
#SBATCH --partition=day
#SBATCH --time=10:00:00
#SBATCH --output=job_output.txt
#SBATCH --cpus-per-task=25

# run optimization
./main

# copy results into results folder
cp nlopt_testfcts.txt ../results/nlopt_testfcts.txt
