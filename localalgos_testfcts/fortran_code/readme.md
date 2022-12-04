# Local algorithms on Test Functions

## Amoeba
comes from original tiktak code, with modification to return {x_min, f_min} even if max number of iterations is reached. 
* Precision: 1.0D-18
* Max iterations: max_evalstep - dim (not sure why again)

## DFNLS
comes from tiktak code (with modification?) [also called boyqa_h]
* precision
    * rhobeg =  (MINVAL(p_bound(:,2)-p_bound(:,1))/2.5_DP)/(4*itratio) with itratio = 1
    * rhoend  = 1.0D-8

* p_ninterppt = 2*dim + 1
* max iterations: max_evalstep
* pnmom = 1

## DFPMIN
comes from Numerical Recipes with some modifications (?)
* Precision: 1.0D-12

## Note:
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