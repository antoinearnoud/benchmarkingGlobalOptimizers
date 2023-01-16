#!/bin/bash
# -----------------------------------------------------------------------------------
#    This script compile and runs the program in tiktak3_incprocess
#     that looks for the global minimum of our income process
#     using Tiktak algorithm.
#
#     It runs a Monte Carlo simulation.
#     The number of implementations of the Monte Carlo simulation (default=10)
#       can be changed below (variable imp).
#
#   This is part of the code of Arnoud, Guvenen, Kleineberg (2019)
# ------------------------------------------------------------------------------------

#Instructions
# ============================================
# in order to minimize the run time, this script creates folders
# and compile the program and runs it in each folder.
# it will create one folder for each algorithm, each maximum number of calls
# to the function specified and for each initial point.

# This script compiles the program, runs it and generates output
# output is stored in monteCarloResultsMin.dat and monteCarloResults.dat files for individual optimization
# script needs to call sbatch_tiktak_a_incomeprocess.sh and sbatch_tiktak_b_incomeprocess.sh

# Parameters & cleaning
#=====================
echo "running Monte Carlo..."

#need to allow modifying txt file, otherwise cannot rewrite
chmod u+x config.txt

#cleanup some stuff
echo "cleanup..."
rm -rf ../results

#BEGIN USER MODIFIABLE DATA
dim=7 #Number of dimensions of the objective; needed to modify config.txt file
#SobolList=(100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000) #List of quantities of sobol points to try; separate elements with space
SobolList=(40 70 100 200 400 800 1600 3000 4000 6000) # for tiktak-d3 and tiktak-nm3
length=${#SobolList[@]}
length=$((length-1)) # length of SobolList (minus 1 because indexation starts at 0)
imp=10   # number of trials per quantity of sobol points
#END USER MODIFIABLE DATA

# create folder results
mkdir ../results

# create and store file with configuration used for benchmarking
echo "number of Sobols (list):" > ../results/monteCarloConfig.dat # create file if file does not exists; replace content if file exists
printf "%s\n" "${SobolList[@]}" >> ../results/monteCarloConfig.dat # print the list of numbers of Sobol points
echo "Local search criterion parameters:" >>  ../results/monteCarloConfig.dat
sed -n '21p' genericParams.f90 >> ../results/monteCarloConfig.dat
sed -n '23p' genericParams.f90 >> ../results/monteCarloConfig.dat

# display info on the terminal
echo "Parameters of the benchmarking:"
cat ../results/monteCarloConfig.dat

# create one folder by algo (DFNLS [bobyqa in the code] and amoeba) and number of Sobol points
for j in $(seq 0 $length); do
  n=${SobolList[$j]}
  mkdir ../results/dfnls_sobol$n
  mkdir ../results/amoeba_sobol$n
done

# create one folder by starting point
for j in $(seq 0 $length); do
  n=${SobolList[$j]}
  for i in $(seq 1 $imp); do
    mkdir ../results/dfnls_sobol$n/startpoint_$i
    mkdir ../results/amoeba_sobol$n/startpoint_$i
  done
done


# copy everything in each folder and change config file appropriately
for i in $(seq 1 $imp); do

  # overwrite starting point on config file
  #(using XX random starting points but always the same point for each step in performance profile)
  for k in $(seq 1 $dim); do
    lineout=$((37+$dim+$k))
    linein=$((($i-1)*$dim+$k))
    sed -i "${lineout}s/.*/`sed -n "${linein}p" init_points.dat`/g" config.txt
  done

  # varying number of Sobol points
  for j in $(seq 0 $length); do
    n=${SobolList[$j]}   # number of sobols generated
    m=$((n/10))	# number of sobols kept

    #################################### dfnls ###########################################
    # copy files in new folder
    cd ..
    cp -r fortran_code results/dfnls_sobol$n/startpoint_$i
    cp -r data results/dfnls_sobol$n/startpoint_$i
    cp -r SWEout results/dfnls_sobol$n/startpoint_$i

    # go to new folder
    cd results/dfnls_sobol$n/startpoint_$i/fortran_code

    # overwrite number of Sobol points generated / kept on config file
    echo "dfnls (bobyqa) dim7 sobol is $n, kept is $m and iter is" $i
    sed -i "28s/.*/7, 297, -1, ${n[@]}, ${m[@]}, 0, -1/g" config.txt

    # go back to original folder
    cd ../../../../fortran_code



    #################################### amoeba ###########################################
    # copy files in new folder
    cd ..
    cp -r fortran_code results/amoeba_sobol$n/startpoint_$i
    cp -r data results/amoeba_sobol$n/startpoint_$i
    cp -r SWEout results/amoeba_sobol$n/startpoint_$i

    # go to new folder
    cd results/amoeba_sobol$n/startpoint_$i/fortran_code

    # overwrite number of Sobol points generated / kept on config file
    echo "amoeba dim7 sobol is $n, kept is $m and iter is" $i
    sed -i "28s/.*/7, 297, -1, ${n[@]}, ${m[@]}, 0, -1/g" config.txt

    #go back to original folder
    cd ../../../../fortran_code


  done
done

cd ../results

# run optimimization in each folder
for i in $(seq 0 $length); do
  for j in $(seq 1 $imp); do
    n=${SobolList[$i]}

    # go to new folder
    cd dfnls_sobol$n/startpoint_$j/fortran_code

    # run minimization with dfnls
    sbatch sbatch_tiktak_b_incomeprocess.sh
    #./GlobalSearch 0 config.txt b &

    # go back to folder results
    cd ../../../

    # go to new folder
    cd amoeba_sobol$n/startpoint_$j/fortran_code

    # run minimization with amoeba. comment this out if only running dfls
    sbatch sbatch_tiktak_a_incomeprocess.sh
    #./GlobalSearch 0 config.txt a &

    # go back to folder results
    cd ../../../

    #wait

  done
done

wait

cd ../fortran_code

echo "done."
echo ""

echo "create file collect_results_from_folders.sh"

# heredoc
# create file with content below and substitute variables (except $(seq ) and loop indexes)
# note: single quote around first occurrence FOE would not substitute variables
# the file is used to collects results from all subfolders once optimization is done
cat << FOE > collect_results_from_folders.sh
  #This file is automatically generated by script_nlopt_incomeprocess.sh
  declare -i n m j imp
  SobolList=(${SobolList[@]})
  for j in \$(seq 0 $length); do
  n=\${SobolList[\$j]}
	touch ../results/min_dfnls_sob\$n.dat
	touch ../results/time_dfnls_sob\$n.dat
	touch ../results/min_amoeba_sob\$n.dat
	touch ../results/time_amoeba_sob\$n.dat
 for i in \$(seq 1 $imp); do
   cat ../results/dfnls_sobol\$n/startpoint_\$i/fortran_code/monteCarloResultsMin.dat >> ../results/min_dfnls_sob\$n.dat
   cat ../results/dfnls_sobol\$n/startpoint_\$i/fortran_code/monteCarloResults.dat >> ../results/time_dfnls_sob\$n.dat
   cat ../results/amoeba_sobol\$n/startpoint_\$i/fortran_code/monteCarloResultsMin.dat >> ../results/min_amoeba_sob\$n.dat
   cat ../results/amoeba_sobol\$n/startpoint_\$i/fortran_code/monteCarloResults.dat >> ../results/time_amoeba_sob\$n.dat
 done
done
FOE

chmod u+x collect_results_from_folders.sh # make file executable

echo ""
echo "all jobs submitted."
echo ""
echo "run collect_results_from_folders.sh once all computations done."
