module OBJECTIVE
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    PRIVATE
    PUBLIC objFun, dfovec, initial0

   interface

      function myFortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: myFortSleep
          integer (c_int), intent (in), VALUE :: seconds
      end function myFortSleep

   end interface

contains

    FUNCTION objFun(theta)
        use genericParams
        implicit none
        REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta
        REAL(DP) :: objFun
        REAL(DP),DIMENSION(p_nmom) :: v_err

        CALL dfovec(p_nx, p_nmom, theta, v_err)
        objFun=v_err(1)


    END FUNCTION objFun

    SUBROUTINE dfovec(n, mv, x, v_err)
        USE nrtype
        USE global
        use genericParams
        IMPLICIT NONE

        INTEGER, INTENT(IN)     :: n, mv
        REAL(DP), DIMENSION(n), INTENT(IN)  :: x
        REAL(DP), DIMENSION(mv),INTENT(OUT) :: v_err
        REAL(DP), DIMENSION(mv) :: val1, val2
        INTEGER :: i
        REAL(DP), DIMENSION(p_nx) :: values

        !Variables for sleep call
        integer(c_int) :: mytime, dur
        
        !Count number of function evaluations
        fe_counter=fe_counter+1
        
        !!Levi13 function; global min is at f(1,1,...,1)=0; LB -10; UB 10
        !v_err(1)=SIN(3.141592653589793D0*x(1)*3.0D0)**2 + &
        !        (SIN(3.141592653589793D0*x(n)*2.0D0)**2+1.0D0)*(x(n)-1.0D0)**2       
        !DO i=1,n-1
        !    v_err(1)=v_err(1)+(SIN(3.141592653589793D0*x(i+1)*3.0D0)**2+1.0D0)*(x(i)-1.0D0)**2
        !ENDDO
        !    v_err(1) = v_err(1) + 1.0D0    ! for relative criterion
        !!end Levi13 function
        
        
        !!griewank function; flobal min is at f(0,0,0,...,0) = 0; LB -100; UB 100
        ! val1=(x(1)**2.0)/200.0
        !do i=2,n
        !val1=val1+(x(i)**2.0)/200.0
        !enddo    
        
        
        !val2=cos(x(1))
        !do i=2,n
        !    val2=val2*cos(x(i)/sqrt(2.0))
        !enddo  
        !v_err=val1-val2+1.0
        !v_err = v_err + 1.0D0 !for relative criterion  
        !!end griewank function
        
        
        !!Rosenbrock function; global min is at f(1,1,...,1)=0; LB -100; UB 100
        !v_err=100*(x(2) - x(1)**2.0)**2.0 + (1.0 - x(1))**2
        !do i=2,n-1
        !    v_err=v_err+(100*(x(i+1) - x(i)**2.0)**2.0 + (1.0 - x(i))**2)
        !enddo 
        !v_err = v_err + 1.0D0 !for relative criterion
        !!end Rosenbrock function
        
        
        !!Rastrigin function; global min is at f(0,0,....,0)=0; LB -5.12; UB 5.12
        v_err=x(1)**2-10.0*cos(2*3.141592653589793D0*x(1))
        do i=2,n
            v_err=v_err+x(i)**2-10.0*cos(2*3.141592653589793D0*x(i))
        enddo
        v_err = v_err + 10.0D0*n    !this is part of the original function
        v_err = v_err + 1.0D0       !for relative criterion
        !!end rastrigin function
        
        
        mytime=1;
        !dur=myFortSleep(mytime);

    END SUBROUTINE dfovec

    SUBROUTINE initial0
        USE global
        IMPLICIT NONE
    END SUBROUTINE initial0
end MODULE objective
