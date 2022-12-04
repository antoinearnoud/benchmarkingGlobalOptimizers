# -----------------------------------------------------------------------------------
#    This script compiles and runs the program in nlopt_testfcts
#     that looks for the global minimum of test functions,
#     using NLopt algorithms:
#           - NLOPT_GN_CRS2_LM,
#           - NLOPT_GN_ISRES,
#           - NLOPT_G_MLSL_LDS,
#           - NLOPT_GN_ESCH,
#           - NLOPT_GD_STOGO
#           - NLOPT_LD_LBFGS
#           - NLOPT_LN_NELDERMEAD
#           - NLOPT_LN_BOBYQA
#     All algorithms (except last three) are followed by one polishing local search (either using Nelder Mead
#       or Bobyqa).
#     This script runs a Monte Carlo simulation.
#     The number of implementations of the Monte Carlo simulation (default=100)
#       can be changed in main.f90 (variable number_iterations).
#
#   This is part of the code for Arnoud, Guvenen, Kleineberg (2020)
# ------------------------------------------------------------------------------------

#Instructions
# ============================================
# in order to minimize the run time, this script can be run with OPENMP.
# You need to link the compiler to the NLOPT library.


# Loading NLOPT library and ifort compiler on Yale cluster
# ============================================
echo "loading NLOPT library and Intel compiler on Grace..."
# ___old Grace style___
#module load Libs/NLopt		# NLOPT with C++ enabled (needfored for STOGO)
#module load Langs/Intel	# ifort compiler
# ___new Grace style___
echo "loading NLopt..."
module load NLopt/2.6.1-GCCcore-7.3.0
echo "   ... done."
echo "loading ifortran..."
module load ifort/2018.3.222-GCC-7.3.0-2.30
echo "   ... done."
echo "libraries loaded"
echo ""

# Deleting old files
# =============================================
echo "deleting old files..."
rm -f *.mod
rm -f *lock
rm -f *~
rm -f -R ../results
mkdir ../results/
echo "   ...old files deleted"
echo ""

# create and store file with configuration used for benchmarking
# echo "number of function evaluations (list):" > ../results/monteCarloConfig.dat # create file if file does not exists; replace content if file exists
# echo "hard coded in main.f90 file" >>  ../results/monteCarloConfig.dat
# echo "Local search criterion parameters:" >>  ../results/monteCarloConfig.dat
# sed -n '42p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '43p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '44p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '45p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '46p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '47p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '48p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '49p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '50p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '51p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '52p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '53p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '54p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '55p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '56p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '57p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '58p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '59p' main.f90 >> ../results/monteCarloConfig.dat
# sed -n '60p' main.f90 >> ../results/monteCarloConfig.dat

# Compile and run
# =============================================
echo "compile files..."
# ___gfortran___
#gfortran -O3 -fopenmp -lnlopt_cxx -lnlopt myparams_counter.f90 main.f90 -o main
gfortran -O3 -fopenmp -lnlopt myparams_counter.f90 testFunctions.f90 main.f90 -o main
# ___ifort___
## new ifort compilers use -fopenmp or -qopenmp (not openmp) i believe.
#ifort -O3 -openmp -lnlopt myparams_counter.f90 main.f90 -o main
echo "   ... compilation done"

echo "submit job"
sbatch sbatch_nlopt_testfcts.sh
echo "   ...job submitted. Results will be directly copied to folder results."
