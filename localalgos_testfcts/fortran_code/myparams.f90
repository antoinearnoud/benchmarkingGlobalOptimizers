MODULE myparams
  ! This module contains variables related to the objective functions

  USE nrtype ! this only useful with income process

  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: dim_list = (/ 2, 10 /)                           ! list of dimensions of the search domain
  INTEGER   :: dim, idim                                                     ! dimension of the domain
  DOUBLE PRECISION, ALLOCATABLE::  lb(:), ub(:)                             ! lower and upper bounds (vectors) of the domain
  DOUBLE PRECISION, ALLOCATABLE:: solution_point(:)                         ! known (analytical) solution point of the minimization
  DOUBLE PRECISION fsol                                                     ! known (analytical) minimum

END MODULE myparams
