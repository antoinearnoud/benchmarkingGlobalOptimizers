# -----------------------------------------------------------------------------------
#    This script tests that the test functions are well written.
#
#   This is part of the code for Arnoud, Guvenen, Kleineberg (2020)
# ------------------------------------------------------------------------------------

# Libraries
#echo "loading NLopt..."
#module load NLopt/2.6.1-GCCcore-7.3.0
#echo "   ... done."

# Compilers
#echo "loading ifortran..."
#module load ifort/2018.3.222-GCC-7.3.0-2.30
#echo "   ... done."

#echo "libraries and compiler loaded"
#echo ""


#gfortran -O3 -fopenmp -lnlopt myparams_counter.f90 testFunctions.f90 test.f90 -o main
echo " start compiling ..."
gfortran -O3 myparams_counter.f90 testFunctions.f90 test.f90 -o main
echo "   ... compilation done"
echo " "

echo "run the tests on test functions ..."
./main
echo " ... end of program."