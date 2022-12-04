module OBJECTIVE
    ! this module contains the objective function (which calls the right test function)
    use nrtype
    use myparams  ! import dim (dimension of the domain space)
    use genericParams, only: p_nx !, p_wspace, p_nmom, p_maxeval, p_iprint, p_maxpoints, p_ninterppt, p_range, p_bound ! added this; dimension of the domain of optimization
    use monteCarloParams ! import funcloop that tells us which test function to run
    use testFunctions

    implicit none
    PRIVATE
    PUBLIC objFun, dfovec, gradf

contains

    FUNCTION objFun(x)
        !use genericParams ! no need: p_nx is directly declare in this file
        implicit none
        REAL(DP),DIMENSION(p_nx),INTENT(IN)     :: x
        REAL(DP)                :: objFun
        REAL(DP)                :: val !,val1, val2, 
        !REAL(DP)              :: pi2
        !INTEGER                 :: i , n

        ! added
        REAL(DP) :: grad(p_nx)
        INTEGER :: need_gradient        

        !n = p_nx  ! used n because that's how the function have been defined in other files
        !pi2 = 3.141592653589793D0

        !Count number of function evaluations
        !counter=counter+1 ! => counter is called in levi13, griewank, rosenbrock, rastrigrin

        need_gradient = 0

        if(funcloop == 1) Then
          call levi13(val, dim, x, grad, need_gradient)
        elseif(funcloop == 2) Then
          call griewank(val, dim, x, grad, need_gradient)
        elseif (funcloop == 3) THEN
          call rosenbrock(val, dim, x, grad, need_gradient)
        elseif (funcloop ==4) Then
          call rastrigin(val, dim, x, grad, need_gradient)
        endif

        objFun = val


    !     if(funcloop == 1) Then                  !!Levi13 function; global min is at f(1,1,...,1)=0; LB -10; UB 10
    !          val = SIN(pi2*x(1)*3.0)**2 + (SIN(pi2*x(n)*2.0)**2+1.0)*(x(n)-1.0)**2   ! this line is different from the one below
    !           do idim=1,n-1
    !              val = val + (SIN(pi2*x(idim+1)*3.0)**2+1.0)*(x(idim)-1.0)**2
    !           enddo
    !    val = val + 1.0 ! add 1 so relative stop criterion on f works

    !    objFun = val

    ! elseif(funcloop == 2) Then          !! griewank
    !    val1 = 0.0
    !    do idim = 1, n
    !       val1 = val1 + (x(idim)**2.0)/200.0
    !    enddo
    !    val2 = 1.0
    !    do idim = 1,n
    !           val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
    !    enddo
    !    val = val1 - val2 + 1.0 ! +1.0 is part of the function
    !    val = val + 1.0         ! add 1 so the relative criterion on f works
    !        objFun=val

    ! elseif(funcloop == 3) then          !! rosenbrock
    !    val = 0.0
    !    do idim = 1,n-1
    !           val=val+(100*(x(idim+1) - x(idim)**2.0)**2.0 + (1.0 - x(idim))**2)
    !    enddo
    !        val = val + 1.0D0 !for relative criterion

    !    objFun = val

    ! elseif(funcloop == 4) then          !! rastrigrin
    !    val = 0.0D0
    !        do idim = 1,n
    !           val = val + x(idim)**2 - 10.0*cos(2*3.141592653589793D0*x(idim))
    !        enddo
    !    val = val + 10.0D0*n !this is part of the original function
    !        val = val + 1.0D0    !for relative criterion

    !    objFun = val

    ! endif

    END FUNCTION objFun

    ! not used in this program
    FUNCTION gradf(x)
        implicit none
        REAL(DP),DIMENSION(p_nx),INTENT(IN)     :: x
          REAL(DP), DIMENSION(p_nx)     :: gradf, grad
        REAL(DP)                :: objFun
        REAL(DP)                :: val1, val2, val
        !REAL(DP)              :: pi2
        !INTEGER                 :: i !, n
        !INTEGER                   :: igrad

        INTEGER :: need_gradient
        !INTEGER :: n ! FIXME
        
        !n = p_nx
        !pi2 = 3.141592653589793D0
        !counter = counter + dim ! => counter is called in levi13, griewank, rosenbrock, rastrigrin

        need_gradient = 1

        if(funcloop == 1) Then
          call levi13(val, dim, x, grad, need_gradient)
        elseif(funcloop == 2) Then
          call griewank(val, dim, x, grad, need_gradient)
        elseif (funcloop == 3) THEN
          call rosenbrock(val, dim, x, grad, need_gradient)
        elseif (funcloop ==4) Then
          call rastrigin(val, dim, x, grad, need_gradient)
        endif
        gradf = grad

    !    if(funcloop == 1) Then                   !!Levi13 function; global min is at f(1,1,...,1)=0; LB -10; UB 10
    !       grad(1) = 3.0*pi2*2.0*sin(3.0*pi2*x(1))*cos(3.0*pi2*x(1)) &
    !         +2.0*(x(1)-1.0)*(1.0+sin(3.0*pi2*x(2))**2.0)
    !       grad(n) = (x(n)-1.0)**2.0*cos(2.0*pi2*x(n))*sin(2.0*pi2*x(n))*2.0*2.0*pi2 &
    !         + 2.0*(x(n)-1.0)*(1+sin(2.0*pi2*x(n))**2.0) &
    !         + ((x(n-1)-1.0)**2.0)*2.0*cos(3.0*pi2*x(n))*sin(3.0*pi2*x(n))*3.0*pi2
    !       if (n>2) then
    !         do igrad = 2, n-1
    !             grad(igrad) = 2.0*(x(igrad)-1.0)*(1.0+sin(3.0*pi2*x(igrad+1))**2.0) &
    !                 + 2.0*cos(3.0*pi2*x(igrad))*sin(3.0*pi2*x(igrad))*3.0*pi2*(x(igrad-1)-1.0)**2.0
    !         enddo
    !       endif

    ! elseif(funcloop == 2) Then          !! griewank

    !        ! i = 1
    !        val2 = 1.0
    !        do idim=2,n
    !          val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
    !        enddo
    !        grad(1) = 2.0*x(1)/200.0 &
    !                 + (1.0/dsqrt(1+0.0D0))*(sin(x(1)/dsqrt(1+0.0D0))) * val2
    !        ! i = n
    !        val2 = 1.0
    !        do idim=1,n-1
    !          val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
    !        enddo
    !        grad(n) = 2.0*x(n)/200.0 &
    !                  + (1.0/dsqrt(n+0.0D0))*(sin(x(n)/dsqrt(n+0.0D0))) * val2

    !        !i in (2,n-1)
    !        if(n>2)then
    !          do igrad=2,n-1
    !             val2 = 1.0
    !             do idim=1,igrad-1
    !               val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
    !             enddo
    !             do idim = igrad+1,n
    !               val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
    !             enddo
    !             grad(igrad) = 2.0*x(igrad)/200.0 &
    !             + (1.0/dsqrt(igrad+0.0D0))*(sin(x(igrad)/dsqrt(igrad+0.0D0))) * val2
    !          enddo
    !        endif


    ! elseif(funcloop == 3) then          !! rosenbrock
    !        grad(1) = -400.0 * x(1) * (x(2) - x(1)**2.0) - 2.0 * (1.0 - x(1))
    !        grad(n) = 200.0 * (x(n) - x(n-1)**2.0)**2.0

    !        if(n>2)then
    !          do igrad=2,n-1
    !            grad(igrad) =  200.0 * (x(igrad) - x(igrad-1)**2.0) &
    !                         - 400.0 * x(igrad) * (x(igrad + 1) -x(igrad)**2.0) &
    !                           - 2.0 * (1.0 - x(igrad))
    !          enddo
    !        endif

    ! elseif(funcloop == 4) then          !! rastrigrin
    !      do igrad = 1, n
    !         grad(igrad) = 2.0 * x(igrad) + 20.0 * 3.141592653589793D0 * sin(2.0 * 3.141592653589793D0 * x(igrad))
    !     enddo
    ! endif
    ! gradf = grad
    END FUNCTION gradf


 SUBROUTINE dfovec(n, mv, x, v_err)
        USE nrtype
        use genericParams ! needded for funcloop
        IMPLICIT NONE

        INTEGER, INTENT(IN)     :: n, mv
        REAL(DP), DIMENSION(n), INTENT(IN)  :: x
        REAL(DP), DIMENSION(mv),INTENT(OUT) :: v_err

        v_err(1) = objFun(x)

    END SUBROUTINE dfovec

end MODULE objective
