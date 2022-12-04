#!/bin/bash
# -----------------------------------------------------------------------------------
#    This script compile and runs the program in tiktak8_testfcts
#     that looks for the global minimum of our test functions
#     using Tiktak algorithm.
#
#   It calls other scripts (script_griewank_2dim, script_griewank_10dim, ...)
#     that are specific to each function and dimension.
#
#   This is part of the code of Arnoud, Guvenen, Kleineberg (2019)
# ------------------------------------------------------------------------------------

# create scripts executables
chmod u+x script_griewank_2dim.sh
chmod u+x script_griewank_10dim.sh
chmod u+x script_rastrigin_2dim.sh
chmod u+x script_rastrigin_10dim.sh
chmod u+x script_levi_2dim.sh
chmod u+x script_levi_10dim.sh
chmod u+x script_rosen_2dim.sh
chmod u+x script_rosen_10dim.sh

# create folders results
rm -rf ../results
mkdir ../results
mkdir ../results/rastrigin2
mkdir ../results/rastrigin10
mkdir ../results/griewank2
mkdir ../results/griewank10
mkdir ../results/levi2
mkdir ../results/levi10
mkdir ../results/rosen2
mkdir ../results/rosen10

# copy code into folders
cd ..
cp -R fortran_code results/rastrigin2
cp -R fortran_code results/rastrigin10
cp -R fortran_code results/griewank2
cp -R fortran_code results/griewank10
cp -R fortran_code results/levi2
cp -R fortran_code results/levi10
cp -R fortran_code results/rosen2
cp -R fortran_code results/rosen10

# run scripts
cd results/rastrigin2/fortran_code
chmod u+x script_rastrigin_2dim.sh
#./script_rastrigin_2dim.sh
sbatch script_rastrigin_2dim.sh

cd ../../..
cd results/rastrigin10/fortran_code
chmod u+x script_rastrigin_10dim.sh
#./script_rastrigin_10dim.sh
sbatch script_rastrigin_10dim.sh

cd ../../..
cd results/griewank2/fortran_code
chmod u+x script_griewank_2dim.sh
#./script_griewank_2dim.sh
sbatch script_griewank_2dim.sh

cd ../../..
cd results/griewank10/fortran_code
chmod u+x script_griewank_10dim.sh
#./script_griewank_10dim.sh
sbatch script_griewank_10dim.sh

cd ../../..
cd results/levi2/fortran_code
chmod u+x script_levi_2dim.sh
#./script_levi_2dim.sh
sbatch script_levi_2dim.sh

cd ../../..
cd results/levi10/fortran_code
chmod u+x script_levi_10dim.sh
#./script_levi_10dim.sh
sbatch script_levi_10dim.sh

cd ../../..
cd results/rosen2/fortran_code
chmod u+x script_rosen_2dim.sh
#./script_rosen_2dim.sh
sbatch ./script_rosen_2dim.sh

cd ../../..
cd results/rosen10/fortran_code
chmod u+x script_rosen_10dim.sh
#./script_rosen_10dim.sh
sbatch script_rosen_10dim.sh
