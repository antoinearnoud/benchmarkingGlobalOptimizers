module myparams_counter
    IMPLICIT NONE
    ! Note: counter needs to be in a module because it is called within the functions which are in testFunctions.f90 module.

    integer :: counter ! counts the number of calls to the function during optimization

    ! declaration used for openmp (need to declare all threadprivate variables as such in their own module)
    !$omp threadprivate(counter)

  end module myparams_counter
