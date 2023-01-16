# Script (bash) to compile globalSearch program (without running it)

echo "remove old files"
rm -f *lock
rm -f *.mod
rm -f monteCarloGOPAmin.dat
rm -f monteCarloGOPA.dat
rm -f resultsMC.dat

echo ""
echo "compile GlobalSearch program"
gfortran nrtype.f90  myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch

echo ""
echo "can run program using ./GlobalSearch 0 config.txt a or ./GlobalSearch 0 config.txt b"
