#!/bin/bash
#SBATCH --job-name=tiktak_tf
#SBATCH --partition=week
#SBATCH --time=06-10:00:00


#echo "rebuild? 1 = yes, 0 = no"
#read rebuild
#if [ $rebuild == 1 ]; then

#generate random starting points using bounds from config file
#gfortran GOPA_init_point.f90 -o initpoint
#./initpoint

gfortran nrtype.f90 stateControl.f90 genericParams.f90 simplex.f90 global.f90 objective_griewank.f90 utilities.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch

#ifort -m64 -g -debug all -heap-arrays nrtype.f90 stateControl.f90 genericParams.f90 simplex.f90 global.f90 objective_griewank.f90 utilities.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch

#fi

if [ $? -ne 0 ]; then
  echo "Errors compiling GlobalSearch"
  exit
fi

rm -f monteCarloResults.dat
rm -f monteCarloResultsMin.dat
rm -f startingPoints.dat


#echo "enter number of iterations and press ENTER"
#iter=1

#echo "enter algortithm (a,d,b) and press ENTER"
#read algo

echo "running tiktak algo..."
echo "     removing lock files"
rm -f *.lock

#need to allow modifying txt file, otherwise cannot rewrite
chmod u+x config_griewank2.txt
chmod u+x config_griewank10.txt

declare -i n m N S NN SS dim lineout linein
#N is max # of steps, S is stepsize, then loop generates n=# Sobols generated (10% of N*S) and then m=#Sobols kept(25% of N*S)

N=10
S=2000
NN=10
SS=8000


#AMOEBA in 10 dimensions varying points generated/kept
dim=10

for i in $(seq 1 100); do
#overwrite starting point (using 100 random starting points but always the same 100 points for each step in performance profile)
	for k in $(seq 1 $dim); do
	lineout=37+$dim+$k
	linein=($i-1)*$dim+$k
	sed -i "${lineout}s/.*/`sed -n "${linein}p" init_points_griewank.dat`/g" config_griewank10.txt
	  done

for j in $(seq 1 10); do
	n=50*$j
	m=$n/10

#echo "amoeba dim10 sobol is $n, kept is $m and iter is" $i
#overwrite number of Sobol points generated / kept
	sed -i "28s/.*/10, 1, -1, ${n[@]}, ${m[@]}, 0, 10/g" config_griewank10.txt
	./GlobalSearch 0 config_griewank10.txt a

done

for j in $(seq 1 20); do
	n=$n+125
	m=$n/10

#echo "amoeba dim10 sobol is $n, kept is $m and iter is" $i
#overwrite number of Sobol points generated / kept
	sed -i "28s/.*/10, 1, -1, ${n[@]}, ${m[@]}, 0, 10/g" config_griewank10.txt
	./GlobalSearch 0 config_griewank10.txt a
done

done


#BOBYQA in 10 dimensions varying points generated/kept
dim=10

for i in $(seq 1 100); do
#overwrite starting point (using 100 random starting points but always the same 100 points for each step in performance profile)
	for k in $(seq 1 $dim); do
	lineout=37+$dim+$k
	linein=($i-1)*$dim+$k
	sed -i "${lineout}s/.*/`sed -n "${linein}p" init_points_griewank.dat`/g" config_griewank10.txt
	  done

for j in $(seq 1 10); do
	n=50*$j
	m=$n/10

#echo "amoeba dim10 sobol is $n, kept is $m and iter is" $i
#overwrite number of Sobol points generated / kept
	sed -i "28s/.*/10, 1, -1, ${n[@]}, ${m[@]}, 0, 10/g" config_griewank10.txt
	./GlobalSearch 0 config_griewank10.txt b

done

for j in $(seq 1 20); do
	n=$n+250
	m=$n/10

#echo "amoeba dim10 sobol is $n, kept is $m and iter is" $i
#overwrite number of Sobol points generated / kept
	sed -i "28s/.*/10, 1, -1, ${n[@]}, ${m[@]}, 0, 10/g" config_griewank10.txt
	./GlobalSearch 0 config_griewank10.txt b
	done

done


# to print the results in a file, use the following after on the previous line:
#> simulations_output_$(date +%Y%m%d).txt

echo "done. check simulations_output file"
