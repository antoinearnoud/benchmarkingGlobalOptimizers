# NLOpt algorithms on Test Functions
Using the NLOpt library to minimize test functions (2 and 10 dimensional domains). 

See NLOpt documentation [here](https://nlopt.readthedocs.io/en/latest/)

## CRES2 (Global + No gradient)
* Name: NLOPT_GN_CRS2_LM
* Precision: 
* Max iterations: 

## ISRES (Global + No Gradient)
* Name: NLOPT_GN_ISRES

## STOGO (Global + Gradient)
* Name: NLOPT_GD_STOGO

## ESCH (Global + No Gradient)
* Name: NLOPT_GN_ESCH

## MLSLS (Global + Local algorithn: NelderMead or BOBYQA)
* Name: NLOPT_G_MLSL_LDS
* Precision local algorithms: either 1.D-3 or 1.D-8 (ftol)

All the above are polished with either Nelder-Mead (from NLopt) or Bobyqa (from NLopt):
* Precision: 1.D-8 (ftol and xtol)

## Local algorithns: LBFGS (similar to dfpmin), NELDERMEAD, BOBYQA
* Names: NLOPT_LN_NELDERMEAD, NLOPT_LN_LBFGS, NLOPT_LN_BOBYQA
* Precision: values by default (FIXME: what are they?)


## Note:
**Counters**
- counter: implemented in each of the test functions, incremented by +1 each time the function is called
- [not here: algo_counter]

**Seed**
For reproducibility the seed is set at iseed = 123456

## RUN THE PROGRAM
Run `bash script_nlopt_testfunctions.sh`: this will run the Monte Carlo implementation and save the results in the (newly created) folder `results`. 

This implementation can be run in parallel (through OpenMP).

It takes approximately XX minutes to run (on Yale cluster) on the login node.
