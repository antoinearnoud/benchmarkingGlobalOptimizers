module monteCarloParams
    !$ use omp_lib
    IMPLICIT NONE
    
    ! Note: counter needs to be in a module because it is called within the functions which are in testFunctions.f90 module.

    integer :: counter            ! counter in the functions (number of calls to functions) -> needs to be in module
    !integer(i4b) :: counter_algo  ! counter in the algo (approximate number of calls)
    !integer :: iloop
    !integer :: funcloop           ! need to be in module -> used in objective for the functions.
    
    !real :: cpu_start, cpu_finish                                             ! to compute time needed for optimization
    !integer, parameter :: number_iterations = 100                             ! number of different starting points for the optimization. used to compute performance statistics
    
    !integer max_evalstep                                                     ! maximum number of function evaluations

    ! declaration used for openmp (need to declare all threadprivate variables as such in their own module)
    !$omp threadprivate(counter)
    
    end module monteCarloParams
    