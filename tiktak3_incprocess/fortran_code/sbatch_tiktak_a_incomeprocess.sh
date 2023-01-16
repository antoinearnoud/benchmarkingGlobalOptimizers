#!/bin/bash
#SBATCH --job-name=tik_inc_a
#SBATCH --partition=week
#SBATCH --time=06-10:00:00

# amoeba (Nelder Mead)
gfortran -O3 nrtype.f90 myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch
./GlobalSearch 0 config.txt a

# SBATCH --output=job_output.txt
