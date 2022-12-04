! -----------------------------------------------------------------------------------
!    This programm looks for the global minimum of test functions
!     using local search algorithms:
!           - Nelder-Mead (called amoeba here),
!           - DFNLS (called -erroneously- bobyqa here),
!           - and dfpmin.
!     It runs a Monte Carlo simulation.
!     The number of implementations of the Monte Carlo simulation (default=100)
!     can be changed with the parameter number_iterations.
!
! ------------------------------------------------------------------------------------

PROGRAM main
  USE nrtype    ! types used in some routines
  USE myparams  ! parameters related to the functions: dimension of the domain, lower and upper bound of the domain, minarg and minimum of the function
  USE minimize  ! contains dfpmin modified from nr; amoeba and dfls from TIKTAK; minor modifications on amoeba (to keep running until max nb of evaluations)
  USE genericParams ! parameters for the solvers
  USE OBJECTIVE, ONLY : objFun, gradf, dfovec
  USE monteCarloParams

  IMPLICIT NONE

  ! Variables Declaration
  !-----------------------
  DOUBLE PRECISION, ALLOCATABLE  :: x(:), x0(:), x_random_mat(:,:)         ! x stores the argmin, x0 is starting point, x_random_mat is matrix of starting points
  DOUBLE PRECISION               :: initpoint                              ! starting point (before it is saved in x_random_mat)
  DOUBLE PRECISION               :: minf                                   ! minf is minimum function value
  INTEGER                        :: algoloop, dim_loop, feloop, jj         ! dummy for loops (algo:algorithms; func:test functions; fe:maximum number of functions evaluations)
  CHARACTER(24)                  :: algorithm , function_name              ! name of algorithm used, function minimized


  
  
  
  
  
   
  ! opens file that will store results
  OPEN(3, file='results_local_algo_all.dat')
  WRITE(3,*) 'thread_no ', 'test_function ', 'algorithm ', 'dimension ', &
       'counter ', 'cpu_time ', 'maxeval ', 'startpt ', 'minf ', 'x1 ', &
       'x2 ', 'x3 ', 'x4 ', 'x5 ', 'x6 ', 'x7 ', 'x8 ', 'x9 ', 'x10 '

  ! looping over dimension parameters (we run the test functions for dimension = 2 and dimension = 10)
  DO dim_loop = 1, 2
     dim = dim_list(dim_loop)
     ALLOCATE(x(dim), x0(dim), x_random_mat(dim, number_iterations))
     ALLOCATE(lb(dim), ub(dim), solution_point(dim))
     ALLOCATE(p_range(dim,2), p_bound(dim,2))

     p_nx = dim          ! needed in ObjFun and other routines (dfnls)

     p_maxpoints = 1 ! 1 is arbitrary.
     p_ninterppt = 2*p_nx+1 !2*nx+1
     p_iprint = 1 !3 ! 0 to 3 for bobyqa_h to be verbose
     p_maxeval = 400*(p_nx+1)
     p_nmom = 1 !number of moments. used for the size of the vector given in dfovec

     ALLOCATE(p_wspace((p_ninterppt+5)*(p_ninterppt+p_nx)+3*p_nx*(p_nx+5)/2)) ! p_wspace is used as working space in bobyqa_h

     PRINT*, "loop over functions starts"
     DO funcloop = 1, 4

        !read in the starting points
        !levi
        IF (funcloop==1) THEN
           OPEN(222, file= "init_points/init_points_levi.dat")
           DO iloop = 1, number_iterations
              DO jj = 1, dim
                 READ(222,*) initpoint
                 x_random_mat(jj,iloop) = initpoint
              ENDDO
           ENDDO
           CLOSE(222)

        !griewank
        ELSE IF (funcloop==2) THEN
           OPEN(221, file= "init_points/init_points_griewank.dat")
           DO iloop = 1, number_iterations
              DO jj = 1, dim
                 READ(221,*) initpoint
                 x_random_mat(jj,iloop) = initpoint
              ENDDO
           ENDDO
           CLOSE(221)

        !rosenbrock
        ELSE IF (funcloop==3) THEN
           OPEN(223, file= "init_points/init_points_rosen.dat")
           DO iloop = 1, number_iterations
              DO jj = 1, dim
                 READ(223,*) initpoint
                 x_random_mat(jj,iloop) = initpoint
              ENDDO
           ENDDO
           CLOSE(223)

        !rastrigin
        ELSE IF (funcloop==4) THEN
           OPEN(224, file= "init_points/init_points_rastrigin.dat")
           DO iloop = 1, number_iterations
              DO jj = 1, dim
                 READ(224,*) initpoint
                 x_random_mat(jj,iloop) = initpoint
              ENDDO
           ENDDO
           CLOSE(224)
        ENDIF

        ! levi13
        IF (funcloop==1) THEN
           function_name = "levi13"
           fsol = 1.0                                  ! because we added 1 to the function
           DO idim=1,dim
              solution_point(idim) = 1.0              ! known solution (used to evaluate success)
              lb(idim) = -10.0D0                      ! lower and upper bounds of the search domain
              ub(idim) = 10.0D0
           ENDDO

           ! griewank
        ELSE IF (funcloop==2) THEN
           function_name = "griewank"
           fsol = 1.0                                  ! because we added 1 to the function
           DO idim=1,dim
              solution_point(idim) = 0.0              ! known solution (used to evaluate success)
              lb(idim) = -100.0                       ! lower and upper bounds of the search domain
              ub(idim) = 95.0                         ! upper bound is 100 but rescale the bounds asymetrically so that it does not directly guess the right point
           ENDDO

           !rosenbrock
        ELSE IF (funcloop==3) THEN
           function_name = "rosenbrock"
           fsol = 1.0
           DO idim=1,dim
              solution_point(idim)=1.0
              lb(idim) = -100.0
              ub(idim) = 100.0
           ENDDO

           !rastrigin
        ELSE IF (funcloop==4) THEN
           function_name = "rastrigin"
           fsol = 1.0
           DO idim=1,dim
              solution_point(idim)=0.0
              lb(idim) = -5.12
              ub(idim) = 5.0 ! theoretically 5.12 but rescale the bounds asymetrically  so that id does not directly guess the right point
           ENDDO
        ENDIF


        ! loop over the different search algorithms
        DO algoloop=1, 3
           IF (algoloop == 1) THEN
              algorithm = 'amoeba'
           ELSEIF (algoloop == 2) THEN
              algorithm = 'dfpmin'
           ELSE
              algorithm = 'dfls'
           ENDIF

           print*, "name of algorithm used: ", algorithm

           ! loop over iterations
           DO iloop = 1, number_iterations

              ! random initial guess point
              x0 = x_random_mat(:,iloop)

              IF (funcloop==1) THEN
                 PRINT*, "function is levi13"
              ELSEIF (funcloop==2) THEN
                 PRINT*, "function is griewank"
              ELSEIF (funcloop==3) THEN
                 PRINT*, "function is rosenbrock"
              ELSEIF (funcloop==4) THEN
                 PRINT*, "function is rastrigrin"
              ENDIF

              ! define maximum number of function evaluations
              max_evalstep=0
              DO feloop = 1, 30
                 x=x0
                 IF (feloop <= 9) THEN
                    max_evalstep =  max_evalstep+50*dim             !100-increments in 2 dim, 500-increments in 10dim
                 ENDIF
                 IF (feloop >= 10) THEN
                    max_evalstep=max_evalstep+dim*10000/20          !1000=increments for 2 dim; 5000 increments for 10 dim
                 ENDIF

                 ! initialize counter of number of calls to the function during search
                 counter = 0       ! counter in the function (increments at every call)
                 counter_algo = 0  ! counter in amoeba

                 ! p_bound and p_range is used to create the first simplex before running amoeaba
                 p_bound(:,1) = lb
                 p_bound(:,2) = ub
                 p_range(:,1) = lb
                 p_range(:,2) = ub

                 CALL CPU_TIME(cpu_start)

                 ! start optimization
                 IF (algoloop == 1) THEN
                  CALL runAmoeba(x, p_tolf_amoeba, max_evalstep-dim)    ! starting point, tolerance, itmax ; p_tolf_amoeba used to be '0.000000000000000001D0'
                    ! [itmax = evlastep - (dim + 1)] because runAmoeba estimates the function once before calling the routine that uses itmax; in Amoeba there will be n+1 (=dim +1) evluation for the simplex
                    minf = objFun(x)    ! evaluation of the function at starting point (informational only)
                    counter = counter-1 ! cancels the call to the function previous line ! FIXME: not sure we want to cancel the last check as the routine does not return the value

                 ELSEIF (algoloop == 2) THEN
                    CALL dfpmin(x, dfpmin_tol,counter_algo,minf,max_evalstep) !0.000000000001D0

                 ELSEIF (algoloop == 3) THEN
                    !search
                    itratio= 1 ! REAL(whichPoint)/REAL(p_maxpoints)
                    IF (itratio>0.5) THEN
                       rhobeg  = (MINVAL(p_bound(:,2)-p_bound(:,1))/2.5_DP)/(4*itratio)
                       rhoend  = 1.0D-8/itratio  !THIS LINE CHANGED BY FATIH, 02/21/2017
                    ELSE
                       rhobeg = MINVAL(p_bound(:,2)-p_bound(:,1))/2.50_DP
                       rhoend  = 1.0D-3
                    END IF
                    CALL bobyqa_h(p_nx,p_ninterppt,x,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,max_evalstep,p_wspace,p_nmom)
                    !call bobyqa_h(p_nx,p_ninterppt,evalParam,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)
                    minf = objFun(x)

                 ENDIF

                 CALL CPU_TIME(cpu_finish)

                 !print results of search
                 !WRITE(*,*) ''
                 !WRITE(*,*) '====RESULTS==='
                 !WRITE(*,*) 'algorithm: ', algorithm
                 !WRITE(*,*) 'function: ', function_name
                 !WRITE(*,*) 'dimension: ', dim
                 !WRITE(*,*) 'iteration: ', iloop
                 !WRITE(*,*) 'maxevalstep', max_evalstep
                 !WRITE(*,*) 'x0 = ', x0
                 !WRITE(*,*) 'found min at ', x
                 !WRITE(*,*) 'min val = ', minf
                 !WRITE(*,*) 'found after (counter)', counter, 'iterations'
                 !WRITE(*,*) 'found in ', cpu_finish - cpu_start, 'cpu seconds'
                 !WRITE(*,*) 'real minval at: ', solution_point
                 !WRITE(*,*) 'real minval = ', fsol
                 !WRITE(*,*) '------------------------------------------------------------------------------'
                 !WRITE(*,*) ''

                 WRITE(3,*) 0, function_name, algorithm, dim, &
                      counter, cpu_finish - cpu_start, max_evalstep, x0(1), minf, x

              END DO ! loop for different maxeval / feloop (run 100 times; here same random draw for each time, before moving to next implementation)
           END DO ! loop for number of implementation (e.g. new random draw)
        ENDDO ! endo algoloop
     ENDDO ! end loop on functions


     DEALLOCATE(lb, ub, solution_point)
     DEALLOCATE(p_range, p_bound)
     DEALLOCATE(x, x0, x_random_mat)
     DEALLOCATE(p_wspace)

  ENDDO ! end loop on dimension
  CLOSE(3)

CONTAINS
  !------------------------------------------------------------
  !               LOCAL ALGORITHMS
  !------------------------------------------------------------
  ! Local algorithms are in minimize.f90
  ! The function runAmoeba below sets the first simplex and call the SUBROUTINE amoeba in minimize.f90

  ! AMOEBA
  SUBROUTINE runAmoeba(amoeba_pt,tolf,itmax)
    !This subroutine executes the amoeba search. It takes the point passed in
    !as a parameter, and generates a simplex centered on it. The points of the
    !simplex are proportional to the number of sobol points generated and the
    !dimension of the parameter space. For example, if we have 2 parameters and
    !100 sobol points, then the simplex distance is approximately 1/10 of the
    !range of each parameter.
    USE simplex, ONLY : simplex_coordinates2
    IMPLICIT NONE
    REAL(DP),     INTENT(in)                 :: tolf
    INTEGER(I4B), INTENT(in)                 :: itmax
    INTEGER(I4B) 				                     :: ia
    REAL(DP), DIMENSION(p_nx), INTENT(inout) :: amoeba_pt
    REAL(DP), DIMENSION(p_nx,p_nx+1)         ::  x
    REAL(DP), DIMENSION(p_nx+1,p_nx)         :: xT
    REAL(DP), DIMENSION(p_nx+1)              :: fnVals
    REAL(DP), DIMENSION(p_nx)                :: temp

    CALL simplex_coordinates2(p_nx, x) ! in simplex.f90; return vector x of coordinates of normalized simplex

    DO ia=1,p_nx+1
       temp = amoeba_pt+x(1:p_nx,ia)*(p_range(:,2)-p_range(:,1))/(DBLE(p_maxpoints)**(1.0_dp/p_nx))
       temp = MAX(MIN(temp,p_bound(:,2)),p_bound(:,1))
       fnVals(ia) = objFun(temp)
       xT(ia,:) = temp
    END DO

    ia=itmax
    CALL amoeba(xT,fnVals,tolf,objFun,counter_algo,itmax)  !in minimize.f90.
                                                          ! Modified by AA-TK to keep running until number of max iterations reached.
    ia = MINLOC(fnVals,1) ! in utilities.f90
    amoeba_pt = xT(ia,:)
  END SUBROUTINE runAmoeba


END PROGRAM main
