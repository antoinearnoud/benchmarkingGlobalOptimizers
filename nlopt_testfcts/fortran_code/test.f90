! -----------------------------------------------------------------------------------
!   Test Program:
!   This programm tests the functions.
!
!   Program written for the paper Benchmarking Global Optimizers, 2020
!    by Antoine Arnoud, Fatih Guvenen and Tatjana Kleineberg
! ------------------------------------------------------------------------------------

program test
  use testFunctions
  implicit none

  integer, parameter :: n = 10
  double precision val, x(n), grad(n), sol
  integer need_gradient, idim, igrad

  need_gradient = 0

  do idim = 1, n
    x(idim) = 1.0
  enddo
  call levi13(val, n, x, grad, need_gradient)
  sol = 1.0
  call assert(val, sol)

  do idim = 1, n
    x(idim) = 0.0
  enddo
  call griewank(val, n, x, grad, need_gradient)
  sol = 1.0
  call assert(val, sol)

  do idim = 1, n
    x(idim) = 1.0
  enddo
  call rosenbrock(val, n, x, grad, need_gradient)
  sol = 1.0
  call assert(val, sol)

  do idim = 1, n
    x(idim) = 0.0
  enddo
  call rastrigin(val, n, x, grad, need_gradient)
  sol = 1.0
  call assert(val, sol)


  CONTAINS
    subroutine assert(a,b)
      implicit none
      double precision :: a
      double precision :: b
      double precision :: sum

      sum = a  - b
      if ( abs (sum) < 1.0D-5 ) then
        write(*,*) 'success'
        !call exit(1)
      else
        write(*,*) 'problem'
        write(*,*) a
        write(*,*) b
      endif
    end subroutine assert

    !include 'testFunctions.f90'

End program test
! END OF FILE
