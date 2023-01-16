# -----------------------------------------------------------------------------------
#    This script compile and runs the program in nlopt_incprocess
#     that looks for the global minimum of our income process
#     using NLopt algorithms:
#           - NLOPT_GN_CRS2_LM (2 variants),
#           - NLOPT_GN_ISRES,
#           - NLOPT_G_MLSL_LDS_3,
#           - NLOPT_GN_ESCH,
#           -
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
# it will create one folder for each maximum number of calls to the function specified
# and for each initial point.
# You need to link the compiler to the NLOPT library.


# Loading NLOPT library and ifort compiler on Yale cluster
# ============================================
echo "loading NLOPT library and Intel compiler on Yale Cluster Grace..."
#module load NLopt/2.6.1-GCCcore-7.3.0	# NLOPT with C++ enabled (needed for STOGO) - loaded in each job run (see sbatch_nlopt_incomeprocess.sh)
#module load ifort/2018.3.222-GCC-7.3.0-2.30	# ifort compiler (we use gfortran, loaded by default)
echo "   ...library loaded"
echo ""

# Deleting old files
# =============================================
echo "deleting old files if necessary..."
rm -f *.mod
rm -f *lock
rm -f monteCarloResultsMin.dat
rm -f monteCarloResults.dat
rm -f startingPoints.dat
rm -f *~
rm -f -R ../results
echo "   ...old files deleted"
echo ""


# creation of folders and compulation and Monte Carlo run
# ==============================================
echo "run Monte Carlo..."

# Careful: this script creates length x imp jobs (see below for values of length and imp), i.e. 350 jobs.
# You need a good scheduler that put some of these jobs on hold in order no to overflow your system.

# variables declaration
declare -i n m j imp
MaxEvalList=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17) # array used to compute the maximum number of function evaluations
length=${#MaxEvalList[@]}
length=$((length-1)) # length of MaxEvalList (minus 1 because indexation starts at 0)
imp=10 #25				# number of implementations (starting points). must be smaller or equal to number_iterations in main.f90

# create folder results
cd ..
mkdir results

# create and store file with configuration used for benchmarking
echo "number of function evaluations (list):" > results/monteCarloConfig.dat # create file if file does not exists; replace content if file exists
printf "%s\n" "${MaxEvalList[@]}" >> results/monteCarloConfig.dat # print the list of numbers of Sobol points
echo "  " >>  results/monteCarloConfig.dat
echo "Local search criterion parameters:" >>  results/monteCarloConfig.dat
echo "----------------------------------" >>  results/monteCarloConfig.dat
sed -n '40p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '41p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '42p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '43p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '44p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '45p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '46p' fortran_code/main.f90 >> results/monteCarloConfig.dat
sed -n '47p' fortran_code/main.f90 >> results/monteCarloConfig.dat

# create results subfolders for each value of maximum number of function evaluations
for j in $(seq 0 $length); do
 n=${MaxEvalList[$j]}
 mkdir results/maxeval$n
done

# copy fortran material from fortran_code to results' subfolders + compile code + run each in each subfolder
# loop over maxeval list
for j in $(seq 0 $length); do
	n=${MaxEvalList[$j]}

  # loop over initial points
	for i in $(seq 1 $imp); do
    # create directory
		mkdir results/maxeval$n/initpoint$i
    # copy fortran code
		cp -r fortran_code results/maxeval$n/initpoint$i
		cp -r SWEout results/maxeval$n/initpoint$i
		cp -r data results/maxeval$n/initpoint$i
		cd results/maxeval$n/initpoint$i/fortran_code
		# change the line in main file: this makes each folder run only its own maxeval and its own initial point
		chmod u+x main.f90
		sed -i "282s/.*/if (.NOT.(feloop == ${n})) cycle/" main.f90
		sed -i "222s/.*/if (.NOT.(iloop == ${i})) cycle/" main.f90
		# compile code (no need of STOGO so library used is -lnlopt not -lnlopt_cxx)
		echo "compiling..."
		gfortran -O3 nrtype.f90 myParams.f90 genericParams.f90 objective.f90 myparamsNLOPT.f90 main.f90 -lnlopt -o main
		echo "done"
		# run optimization
		echo "submit job"
		sbatch sbatch_nlopt_incomeprocess.sh 		# syntax for IBM clusters (old Yale clusters): bsub -q week -W 150:00 ./main &
		echo "job submitted"
		# get out
		cd ../../../..
	done
done

echo ""
echo "done."
echo ""

echo "create file collect_results_from_folders.sh"
cd fortran_code

# heredoc
# create file with content below and substitute variables (except $(seq ))
# note: single quote around first occurance FOE would not substitute variables
# the file is used to collects results from all subfolders once optimization is done
cat << FOE > collect_results_from_folders.sh
#This file is automatically generated by script_nlopt_incomeprocess.sh
#It should be run once all the jobs for function otpimization are done running.

  declare -i n m j imp
  MaxEvalList=(${MaxEvalList[@]})
  touch ../results/nlopt_incprocess.dat
  for j in \$(seq 0 $length); do
    n=\${MaxEvalList[\$j]}
    for i in \$(seq 1 $imp); do
      cat ../results/maxeval\$n/initpoint\$i/fortran_code/results_nlopt_incomeprocess.txt >> ../results/nlopt_incprocess.dat
    done
  done
FOE

chmod u+x collect_results_from_folders.sh # make file executable -> can be run simply by ./collect_results_from_folders.sh

echo ""
echo "all jobs submitted."
echo ""
echo "run collect_results_from_folders.sh once all computations done."
