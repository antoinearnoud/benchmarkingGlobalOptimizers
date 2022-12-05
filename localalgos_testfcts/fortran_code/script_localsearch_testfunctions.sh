# -----------------------------------------------------------------------------------
#    This script compile and runs the program in localalgos_testfcts
#     that looks for the global minimum of test functions
#     using local search algorithms:
#           - Nelder-Mead (also called amoeba),
#           - DFNLS (called bobyqa in the main program),
#           - and dfpmin.
#     It runs a Monte Carlo simulation.
#     The number of implementations of the Monte Carlo simulation (default=100)
#       can be changed in the main.f90 file (variable number_iterations).
#
#   This is part of the code of Arnoud, Guvenen, Kleineberg (2019)
# ------------------------------------------------------------------------------------

#Instructions
# ============================================
# No extra library required. This script is not using openmp.

# Deleting old files
# =============================================
echo "delete old files"
rm -f main
rm -f *.mod
rm -f *~
rm -rf ../results
mkdir ../results
echo ""

# Compilation (gfortran)
# ==============================================
echo "compiling..."
# using gfortran
gfortran -O3  nrtype.f90 nrutil.f90 monteCarloParams.f90 genericParams.f90 simplex.f90 utilities.f90 myparams.f90 testFunctions.f90 objective.f90 minimize.f90 main.f90  -o main
# using ifort
#ifort -O3 nrtype.f90 nrutil.f90 monteCarloParams.f90 genericParams.f90 simplex.f90 utilities.f90 myparams.f90 testFunctions.f90 objective.f90 minimize.f90 main.f90 -o main
echo ""
echo "...compilation done"
echo ""

# Minimization
# ==============================================
echo "running program: minimizing test functions with local algorithms..."
sleep 2
#./main
#echo ""
#cp results_local_algo_all.dat ../results/results_local_algo_all.dat
#echo "Local aglo optimization on test functions done. See results_local_algo_all.txt file in results directory."
