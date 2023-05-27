module monteCarloParams
! This module contains counters, the varaible storing the maximum number of evaluations, and the number of iterations. 

USE NRTYPE
IMPLICIT NONE

    integer :: counter            ! counter in the functions (number of calls to functions)
    integer(i4b) :: counter_algo  ! counter in the algo (approximate number of calls)
    integer :: iloop               ! loop over iteration
    integer :: funcloop            ! loop over the functions used; needs to be global var because used in objective.f90

    real :: cpu_start, cpu_finish                                             ! to compute time needed for optimization
    integer, parameter :: number_iterations = 100                             ! number of different starting points for the optimization. used to compute performance statistics

    integer max_evalstep                                                     ! maximum number of function evaluations

end module monteCarloParams
