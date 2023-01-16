#!/bin/bash
#SBATCH --job-name=nlopt_incproc
#SBATCH --partition=week
#SBATCH --time=06-10:00:00

# This is the job that will be run.
module load NLopt/2.6.1-GCCcore-7.3.0
./main
