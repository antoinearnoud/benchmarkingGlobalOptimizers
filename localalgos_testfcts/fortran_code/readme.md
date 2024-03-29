# Local algorithms on Test Functions

This folder contains the code to run three local algorithms (Nelder-MEad, DFNLS, DFPMIN) on 4 test functions (Levi, Griewank, Rosenbrock, Rastrigin).

The files are the following:
- `testFunctions.f90` contains the code of each of the test functions
- `objective.f90` is a general function which calls the right test function depending on the loop index
- `myparams.f90` contains parameters related to the objective function (dimensions, theoretical minimum, theoretical argmin, lower and upper bounds of the domain)
- `minimize.f90` contains the code of the minimization routines: amoeba (Nelder-Mead), bobyqa_h (DFNLS), and DFPMIN (see below)
- `genericParams.f90` contains variables needed for the optimization (tolerance, range of optimization domain, size of the simplex for Amoeba, etc.)
- `monteCarloParams.f90` contains dummies (for loops in the main program), maximum number of evaluations, and the number of iterations of the Monte Carlo experiment
- `simplex.f90`, `nrtype.f90`, `nrutil.f90`, `utilities` are auxiliary functions (create a simplex, define types, etc.)
- `main.f90` is the main program that runs the opimizations and saves the results

## Amoeba
the code comes from the original Tiktak code, with modification to return {x_min, f_min} even if max number of iterations is reached. 
* Precision: 1.0D-18
* Max iterations: max_evalstep - dim

## DFNLS
the code comes from original Tiktak code (with modifications) (the optimization is also called boyqa_h)
* precision
    * `rhobeg =  (MINVAL(p_bound(:,2) - p_bound(:,1)) / 2.5_DP) / (4*itratio)` with `itratio = 1`
    * `rhoend  = 1.0D-8`

* `p_ninterppt = 2*dim + 1`
* max iterations: max_evalstep
* pnmom = 1

## DFPMIN
comes from Numerical Recipes with some modifications
* Precision: 1.0D-12

## Notes:
**Counters**
- counter: implemented in each of the test functions, incremented by +1 each time the function is called
- algo_counter: argument in amoeba minimization and dfpmin minimization; 
    * used in amoeba to stop once maximum number of iterations is reached ; 
    * used in dfpmin to stop once maximum number of iterations is reached; 
    * note that in DFNLS (bobyqa_h) the feature is directly implemented and the algorithm does not go over number of max iterations; 
    * not reported in the result file

## RUN THE PROGRAM
Run `bash script_localsearch_testfunctions.sh`: this will run the Monte Carlo implementation and save the results in the (newly created) folder `results`. 

This implementation does not run in parallel (neither through OpenMP nor through dividing parameters into different groups and calling different jobs).

It takes approximately 10 minutes to run (on Yale cluster) on the login node.
