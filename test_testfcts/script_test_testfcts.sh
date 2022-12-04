# -----------------------------------------------------------------------------------
#    This script tests that the test functions are well written.
#
#   This is part of the code for Arnoud, Guvenen, Kleineberg (2020)
# ------------------------------------------------------------------------------------


echo " start compiling ..."
gfortran -O3 myparams_counter.f90 testFunctions.f90 test.f90 -o main
echo "   ... compilation done"
echo " "

echo "run the tests on test functions ..."
./main
echo " ... end of program."
