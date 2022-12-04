! This file contains the test functions as subroutines.
! Functions: levi13, griewank, rosenbrock, rastrigin
! Note: the file imports the variable "counter" from genericParams to report number of calls

MODULE testFunctions
  IMPLICIT NONE

  Contains

   ! Each subroutine takes the following arguments:
   !  val is output value, 
   !  n is dimension of domain, 
   !  x is where f is evaluated, 
   !  grad is the gradient,
   !  need_gradient = 0 or 1 if one wants to compute the gradient

  subroutine levi13(val, n, x, grad, need_gradient)  
        use genericParams, only: counter
        double precision val, x(n), grad(n)
        integer n, need_gradient, idim, igrad
        double precision pi2

        pi2 = 3.141592653589793D0

        counter = counter + 1

        if (need_gradient.ne.0) then

            counter = counter + 1

            grad(1) = 3.0*pi2*2.0*sin(3.0*pi2*x(1))*cos(3.0*pi2*x(1)) &
                +2.0*(x(1)-1.0)*(1.0+sin(3.0*pi2*x(2))**2.0)
            grad(n) = (x(n)-1.0)**2.0*cos(2.0*pi2*x(n))*sin(2.0*pi2*x(n))*2.0*2.0*pi2 &
                + 2.0*(x(n)-1.0)*(1+sin(2.0*pi2*x(n))**2.0) &
                + ((x(n-1)-1.0)**2.0)*2.0*cos(3.0*pi2*x(n))*sin(3.0*pi2*x(n))*3.0*pi2

            if (n>2) then
                do igrad = 2, n-1
                    grad(igrad) = 2.0*(x(igrad)-1.0)*(1.0+sin(3.0*pi2*x(igrad+1))**2.0) &
                        + 2.0*cos(3.0*pi2*x(igrad))*sin(3.0*pi2*x(igrad))*3.0*pi2*(x(igrad-1)-1.0)**2.0
                enddo
            endif

        endif

        val = SIN(pi2*x(1)*3.0)**2 + (SIN(pi2*x(n)*2.0)**2+1.0)*(x(n)-1.0)**2   ! this line is different from the one below
        do idim=1,n-1
                val = val + (SIN(pi2*x(idim+1)*3.0)**2+1.0)*(x(idim)-1.0)**2
        enddo
        val = val + 1.0 ! add 1 so relative stop criterion on f works
    end subroutine levi13


    subroutine griewank(val, n, x, grad, need_gradient)
           use genericParams, only: counter
           double precision val, x(n), grad(n), val1, val2
           integer n, need_gradient, idim, igrad

           counter = counter + 1

           if (need_gradient.ne.0) then
               counter = counter + 1

               ! i = 1
               val2 = 1.0
               do idim=2,n
                 val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
               enddo
               grad(1) = 2.0*x(1)/200.0 &
                        + (1.0/dsqrt(1+0.0D0))*(sin(x(1)/dsqrt(1+0.0D0))) * val2
               ! i = n
               val2 = 1.0
               do idim=1,n-1
                 val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
               enddo
               grad(n) = 2.0*x(n)/200.0 &
                         + (1.0/dsqrt(n+0.0D0))*(sin(x(n)/dsqrt(n+0.0D0))) * val2

               !i in (2,n-1)
               if(n>2)then
                 do igrad=2,n-1
                    val2 = 1.0
                    do idim=1,igrad-1
                      val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
                    enddo
                    do idim = igrad+1,n
                      val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
                    enddo
                    grad(igrad) = 2.0*x(igrad)/200.0 &
                    + (1.0/dsqrt(igrad+0.0D0))*(sin(x(igrad)/dsqrt(igrad+0.0D0))) * val2
                 enddo
               endif

           endif

           val1 = 0.0
           do idim = 1, n
               val1 = val1 + (x(idim)**2.0)/200.0
           enddo
           val2 = 1.0
           do idim = 1,n
                   val2 = val2 * cos(x(idim)/sqrt(idim+0.0D0))
           enddo
           val = val1 - val2 + 1.0 ! +1.0 is part of the function
           val = val + 1.0         ! add 1 so the relative criterion on f works

    end subroutine griewank


    subroutine rosenbrock(val, n, x, grad, need_gradient)
           use genericParams, only: counter
           double precision val, x(n), grad(n)
           integer n, need_gradient, idim, igrad
           counter = counter + 1

           if (need_gradient.ne.0) then
               counter = counter + 1

               grad(1) = -400.0 * x(1) * (x(2) - x(1)**2.0) - 2.0 * (1.0 - x(1))
               grad(n) = 200.0 * (x(n) - x(n-1)**2.0)**2.0

               if(n>2)then
                 do igrad=2,n-1
                   grad(igrad) =  200.0 * (x(igrad) - x(igrad-1)**2.0) &
                                - 400.0 * x(igrad) * (x(igrad + 1) -x(igrad)**2.0) &
                                - 2.0 * (1.0 - x(igrad))
                 enddo
               endif
           endif

           val = 0.0
           do idim = 1,n-1
               val=val+(100*(x(idim+1) - x(idim)**2.0)**2.0 + (1.0 - x(idim))**2.0)
           enddo
           val = val + 1.0D0 !for relative criterion

    end subroutine rosenbrock


    subroutine rastrigin(val, n, x, grad, need_gradient)
        use genericParams, only: counter
        double precision val, x(n), grad(n)
        integer n, need_gradient, idim, igrad
        counter = counter + 1

        if (need_gradient.ne.0) then
            counter = counter + 1
            do igrad = 1, n
                grad(igrad) = 2.0 * x(igrad) + 20.0 * 3.141592653589793D0 * sin(2.0 * 3.141592653589793D0 * x(igrad))
            enddo
        endif

        val = 0.0D0
            do idim = 1,n
                val = val + x(idim)**2 - 10.0*cos(2*3.141592653589793D0*x(idim))
            enddo
        val = val + 10.0D0*n !this is part of the original function
        val = val + 1.0D0    !for relative criterion
    end subroutine rastrigin

END MODULE
