! -----------------------------------------------------------------------------------
!    This programm looks for the global minimum of our income process
!     using NLopt algorithms.
!     It runs a Monte Carlo simulation.
!     The number of implementations (default=10) can be changed with the parameter
!     number_iterations.
!
!    Program written for the paper Benchmarking Global Optimizers, 2019
!     by Antoine Arnoud, Fatih Guvenen and Tatjana Kleineberg
! ------------------------------------------------------------------------------------

! Make sure that line 226 has a comment indicating it is line 226. If not, need to change the script script_nlopt_incomeprocess.sh
! Similarly for line 181.
PROGRAM main
  !$ use omp_lib          ! omp_lib is used for parallelization with openmp
    USE nrtype
    USE genericParams       ! used in objFun (through p_nx - see objective.f90)
    USE myParams
    USE OBJECTIVE, only : initial0
    USE myparamsNLOPT, only: counter

    IMPLICIT NONE

    external incomeprocess					       ! income process function (defined below)
    integer*8 opt, polish_opt, local_opt
    double precision, allocatable:: x(:), x0(:), polish_x(:), x_random(:), x_random_mat(:,:)
    double precision :: minf
    integer ires
    integer iloop, algoloop, funcloop, feloop, polishloop
    integer idim
    integer converged, max_evalstep
    double precision val, initpoint
    double precision, allocatable:: grad(:)
    character(35) algorithm, algorithm_polish, function_name
    integer       thread_no, nthreads          ! used with openmp
    double precision stopval
    integer jj

    real :: cpu_start, cpu_finish                     ! computes time needed for optimization
    integer, parameter :: number_iterations = 10      !number of implementations (Monte Carlo). used to compute performance profiles
    integer, parameter :: max_eval_mlsl_local = 5000  ! maximum number of evaluations (stopping criterion).
    double precision :: ftol_mlsl_local
    double precision, parameter :: ftol_mlsl_local_3 = 1.D-3
    double precision, parameter :: ftol_mlsl_local_8 = 1.D-8
    integer :: maxeval_polish ! for example max_evalstep / 5 or max_evalstep/2
    double precision, parameter :: ftol_polish = 1.D-8
    double precision, parameter :: xtol_polish = 1.D-8
    integer   :: dim                                ! dimension of the problem
    double precision, allocatable::  lb(:), ub(:)   ! lower and upper bounds (vectors) of the domain
    double precision, allocatable:: solution_point(:) ! known solution of the minimization
   integer, parameter               :: iseed = 123456




    include 'nlopt.f'

    print*,"calling initial0 (load data for income process)"
    call initial0

    p_nx = 7
    p_nmom = 297
    p_maxeval = -1
    p_ninterppt = 2*p_nx+1
    p_qr_ndraw = 50
    p_maxpoints = 4
    p_searchType = 0
    p_lotteryPoints = -1

    IF (allocated(p_range) .eqv. .FALSE.) THEN
        allocate(p_range(p_nx,2))
    END IF
    IF (allocated(p_init) .eqv. .FALSE.) THEN
        allocate(p_init(p_nx))
    END IF
    IF (allocated(p_bound) .eqv. .FALSE.) THEN
        allocate(p_bound(p_nx,2))
    END IF






    open(3, file='results_nlopt_incomeprocess.txt')
    write(3,*)    'thread_no ', 'number_of_threads ', 'seed ', 'test_function ', 'dimension ', &
                  'algorithm ', 'algorithm_polish ',  &
                  'maxeval ', 'maxeval_mlsl_local ',  'maxeval_polish ', &
                  'ftol_mlsl_local ', 'ftol_polish ', 'xtol_polish ' , 'startpt ', &
                  'number_iterations ', 'iloop ', 'feloop ', &
                  'converged ', 'ires ','cpu_time ', 'counter ',  &
                  'minf ',  &
                  'x1 ', 'x2 ', 'x3 ', 'x4 ', 'x5 ', 'x6 ', 'x7 ' !, 'x8 ', 'x9 ', 'x10 '

  ! saving results beefore polishing phase (in case need to re-run only the polishing phase)
  open(4, file='nlopt_incomeprocess_pre_polish.txt')
  write(4,*)    'thread_no ', 'number_of_threads ', 'seed ', 'test_function ', 'dimension ', &
                'algorithm ', 'algorithm_polish ',  &
                'maxeval ', 'maxeval_mlsl_local ',  'maxeval_polish ', &
                'ftol_mlsl_local ', 'ftol_polish ', 'xtol_polish ' , 'startpt ', &
                'number_iterations ', 'iloop ', 'feloop ', &
                'converged ', 'ires ', 'cpu_time ', 'counter ',  &
                'minf ',  &
                'x_prepolish_1 ', 'x_prepolish_2 ', 'x_prepolish_3 ', 'x_prepolish_4 ', 'x_prepolish_5 ', &
                'x_prepolish_6 ', 'x_prepolish_7 ' !, 'x_prepolish_8 ', 'x_prepolish_9 ', 'x_prepolish_10 '

    dim = 7
    function_name = 'income_process'
    allocate(lb(dim), ub(dim), solution_point(dim))
    allocate(x(dim), polish_x(dim), x0(dim), x_random_mat(dim, number_iterations))
    allocate(grad(dim))

    open(222, file= "init_points_incomeprocess.dat")
    do iloop = 1, number_iterations

	call nlosr(iseed)

        do jj = 1, dim
             read(222,*) initpoint
             x_random_mat(jj,iloop) = initpoint
        enddo
    enddo
    close(222)

    lb(1) = 0.05                        ! lower and upper bounds of the search domain
    lb(2) = 0.6
    lb(3) = .1
    lb(4) = -1.
    lb(5) = 0.01
    lb(6) = 0.01
    lb(7) = 0.01

    ub(1) = .5
    ub(2) = .9
    ub(3) = 1.
    ub(4) = -.1
    ub(5) = 1.
    ub(6) = 1.
    ub(7) = 1.

    p_range(:,1) = lb
    p_range(:,2) = ub


   lb(1) = 0.001                        ! lower and upper bounds of the search domain
   lb(2) = 0.1
   lb(3) = .1
   lb(4) = -1.
   lb(5) = 0.001
   lb(6) = 0.001
   lb(7) = 0.001

   ub(1) = 1.
   ub(2) = .99
   ub(3) = 1.
   ub(4) = -.1
   ub(5) = 1.
   ub(6) = 1.
   ub(7) = 1.

   p_bound(:,1) = lb
   p_bound(:,2) = ub

   DO jj=1,p_nx
       p_init(jj) = 0.0
   END DO

        ! loop over the different search algorithms (not using algoloop = 6 and algoloop = 7)
        do algoloop=1,8 !9,11 !5
            !if (algoloop .neq. 1) then cycle    ! to keep on line 148 or change script
	          if(algoloop == 4) then ! skip STOGO (not working with income process)
		           cycle
	          endif
            if (algoloop == 1) then
                algorithm = 'NLOPT_GN_CRS2_LM'
            else if (algoloop == 2) then
                algorithm = 'NLOPT_GN_ISRES'
            else if (algoloop == 3) then
               algorithm = 'NLOPT_G_MLSL_LDS_NELDERMEAD_3'
            else if (algoloop == 4) then
               algorithm = 'NLOPT_GD_STOGO'
            else if (algoloop == 5) then
               algorithm = 'NLOPT_GN_ESCH'
            else if (algoloop == 6) then
               algorithm = 'NLOPT_G_MLSL_LDS_NELDERMEAD_8' !'NLOPT_LN_NELDERMEAD'    !amoeba
            else if (algoloop == 7) then
               algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_3'         !similar to dfpmin
             else if (algoloop == 8) then
               algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_8'
             else if (algoloop == 9) then
               algorithm = 'NLOPT_LD_LBFGS'              !similar to dfpmin
             else if (algoloop == 10) then
               algorithm = 'NLOPT_LN_NELDERMEAD'
             else if (algoloop == 11) then
               algorithm = 'NLOPT_LN_BOBYQA'
            endif

print*, "got here with algo:"
write(*,*) algorithm

            ftol_mlsl_local = 0.0 ! set value 0.0 for algos other than MLSL

            ! loop over iterations (each starts with different initial point)
            !$omp parallel do if(number_iterations.gt.1) default(none), &
            !$OMP& private(max_evalstep, x0, x, polish_x, minf, ires, converged, feloop), &
            !$OMP& private(opt, polish_opt, local_opt, thread_no), &
            !$OMP& shared(x_random_mat, algoloop, dim, lb, ub, funcloop), &
            !$OMP& private(maxeval_polish), &
            !$OMP& shared(function_name, algorithm, solution_point, nthreads), &
            !$OMP& shared(p_nx, p_nmom, p_qr_ndraw, p_ninterppt, p_init, p_bound, p_range ), &
            !$OMP& shared(p_maxeval, p_maxpoints, p_lotteryPoints)
            do iloop=1, number_iterations
                nthreads = -1 ! USED WHEN NOT USING OMP; overwritten by line below when using OMP
                !$ nthreads = OMP_GET_NUM_THREADS()
                thread_no = 0 ! if not using openmp, initialize the thread number
                !$ thread_no = OMP_GET_THREAD_NUM()
                !$ print*, thread_no

print*, "got here with iteration:"
write(*,*) iloop

                !if (iloop .neq. 1) then cycle ! to keep on line 181 or change script. the script creates one different folder for each iteration
                call cpu_time(cpu_start)
                converged = 0               ! will be 1 when the algorithm has converged (not the case for example when maximum number of evaluations is reached)
                ! random initial guess point (only difference between iterations)
                x0 = x_random_mat(:,iloop)
                ! creates the opt object (indicates the algorithm used and the dimensionality)
                if (algoloop==1) then
                    call nlo_create(opt, NLOPT_GN_CRS2_LM, dim)
                else if (algoloop==2) then
                    call nlo_create(opt, NLOPT_GN_ISRES, dim)
                else if (algoloop==3) then
                    call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
                    call nlo_create(local_opt, NLOPT_LN_NELDERMEAD, dim)
                    ftol_mlsl_local = ftol_mlsl_local_3
                    call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local)  ! it is ok to start with a high tolerance for local search and do a final local search at the end... (they say)
                    !call nlo_set_xtol_rel(ires, local_opt, 1.D-3) ! MLSL defaults to ftol_rel=10−15 and xtol_rel=10−7 for the local searches
                    call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
                    call nlo_set_local_optimizer(ires, opt, local_opt)
                else if (algoloop==4) then
                    call nlo_create(opt, NLOPT_GD_STOGO, dim)
                else if (algoloop==5) then
                    call nlo_create(opt, NLOPT_GN_ESCH, dim)
                else if (algoloop==6) then
                  call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
                  call nlo_create(local_opt, NLOPT_LN_NELDERMEAD, dim)
                  ftol_mlsl_local = ftol_mlsl_local_8
                  call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_8)  ! it is ok to start with a high tolerance for local search and do a final local search at the end... (they say)
                  call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
                  call nlo_set_local_optimizer(ires, opt, local_opt)
                  !call nlo_create(opt, NLOPT_LN_NELDERMEAD, dim)
                else if (algoloop==7) then
                  call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
                  call nlo_create(local_opt, NLOPT_LN_BOBYQA, dim)
                  ftol_mlsl_local = ftol_mlsl_local_3
                  call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_3)  ! it is ok to start with a high tolerance for local search and do a final local search at the end... (they say)
                  !call nlo_set_xtol_rel(ires, local_opt, 1.D-3) ! MLSL defaults to ftol_rel=10−15 and xtol_rel=10−7 for the local searches
                  call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
                  call nlo_set_local_optimizer(ires, opt, local_opt)
                else if (algoloop==8) then
                  call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
                  call nlo_create(local_opt, NLOPT_LN_BOBYQA, dim)
                  ftol_mlsl_local = ftol_mlsl_local_8
                  call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_8)  ! it is ok to start with a high tolerance for local search and do a final local search at the end... (they say)
                  call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
                  call nlo_set_local_optimizer(ires, opt, local_opt)
                end if
                ! set lower and upper bounds
                call nlo_set_lower_bounds(ires, opt, lb)
                call nlo_set_upper_bounds(ires, opt, ub)
                ! set the objective function to minimize (could call nlo_set_max_objective)
                call nlo_set_min_objective(ires, opt, incomeprocess, 0)
                max_evalstep=0
                do feloop = 1, 17        ! 1000 2000 3000 4000 5000 6000 7000 8000 9000 10,000 and then 20,000 30,000 40,000 50,000 60,000 70,000 80,000
                    x=x0                      !reset to starting point: x is now the output from last optimization
                    if (feloop <= 10) then
                        max_evalstep=max_evalstep + 1000
                      endif
                    if (feloop >= 11) then
                        max_evalstep= max_evalstep +  10000          !1000=increments; later do 50x1000 increments... should probably differ across dim...25,000 and 100,000?
                    endif
                    !if (feloop .neq. 1) then cycle   ! to keep on line 226 or change script

            print*, "maxevalstep :"
            write(*,*) max_evalstep
                    call nlo_set_maxeval(ires, opt, max_evalstep)

    !define other stopping criteria to avoid early-stops of local algos. better results without these criteria, which are de-activated by default.
        !call nlo_set_ftol_rel(ires, opt, 1.D-20)
        !call nlo_set_xtol_rel(ires, opt, 1.D-20)
        !call nlo_set_ftol_abs(ires, opt, 1.D-20)
        !call nlo_set_xtol_abs(ires, opt, 1.D-20)

    ! initialize counter of number of calls to the function during search
        counter = 0

  ! start optimization
print*, "starting optimization"

        call nlo_optimize(ires, opt, x, minf) ! optimizes the object opt. returns x and minf. ires is message or error.
print*, "done main optimization"

        call cpu_time(cpu_finish)

                !print results of search
                if (ires.lt.0) then
                    write(*,*) 'nlopt failed! error code: ', ires

                    converged = -1
                    algorithm_polish = 'NONE'
                    maxeval_polish = -999999
                    minf = -999999

                polish_x = (/ -999999, -999999, -999999, -999999, -999999, -999999, -999999 /)

                x = polish_x
                ! write results in file (important information is: converged and ires)
                write(3,*)  thread_no, nthreads, iseed, function_name, dim, &
                            algorithm, algorithm_polish, &
                            max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                            ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                            number_iterations, iloop, feloop, &
                            converged, ires, cpu_finish - cpu_start, counter, &
                            minf, &
                            polish_x
                ! write interemediate result (before polishing part - only difference with above is last variable x vs polish_x)
                write(4,*)  thread_no, nthreads, iseed, function_name, dim, &
                            algorithm, algorithm_polish, &
                            max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                            ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                            number_iterations, iloop, feloop, &
                            converged, ires, cpu_finish - cpu_start, counter, &
                            minf, &
                            x

                    ! write(3,*)  thread_no, iloop, feloop, function_name, nsim, nhhsim, algorithm, dim, &
                    !                      converged, counter, cpu_finish - cpu_start, &
                    !                      max_evalstep, x0(1), ' 999999', ' 999999', ' 999999'
                else
                    write(*,*) 'algorithm: ', algorithm
                     write(*,*) ' function: ', 'income process'
                    write(*,*) 'iteration: ', iloop
!                    write(*,*) 'x0 = ', x0
                    write(*,*) 'found min at ', x
                    write(*,*) 'min val = ', minf
                    write(*,*) 'found after ', counter, 'iterations'
                    write(*,*) 'found in ', cpu_finish - cpu_start, 'cpu seconds'
                    write(*,*) 'ires code: ', ires
                    write(*,*) '------------------------------------------------------------------------------'
                    write(*,*) ''


                maxeval_polish = max_evalstep/5

                if (algoloop .LE. 8) then
                  ! loop over two types of polishing (NLOPT-BOBYQA and NLOPT-NELDER_MEAD)
                  do polishloop = 1, 2
                    polish_x = x                   !start at best-point returned so far from the global algo
                    if (polishloop == 1) then
                      algorithm_polish = 'NELDER_MEAD_POLISH'
                      call nlo_create(polish_opt, NLOPT_LN_NELDERMEAD, dim)
                    else
                      algorithm_polish = 'BOBYQA_POLISH'
                      call nlo_create(polish_opt, NLOPT_LN_BOBYQA, dim)
                    endif

                    call nlo_set_lower_bounds(ires, polish_opt, lb)
                    call nlo_set_upper_bounds(ires, polish_opt, ub)
                    call nlo_set_ftol_rel(ires, polish_opt, ftol_polish)
                    call nlo_set_xtol_rel(ires, polish_opt, xtol_polish)
                    call nlo_set_maxeval(ires, polish_opt, maxeval_polish)

                    call nlo_set_min_objective(ires, polish_opt, incomeprocess, 0)
                    call nlo_optimize(ires, polish_opt, polish_x, minf)

                    call cpu_time(cpu_finish)

                    converged = 1

                    ! write results in file
                    write(3,*)  thread_no, nthreads, iseed, function_name, dim, &
                                algorithm, algorithm_polish, &
                                max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                                ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                                number_iterations, iloop, feloop, &
                                converged, ires, cpu_finish - cpu_start, counter, &
                                minf,   &
                                polish_x

                  ! write interemediate result (before polishing part - only difference with above is last variable x vs polish_x)
                  write(4,*)  thread_no, nthreads, iseed, function_name, dim, &
                              algorithm, algorithm_polish, &
                              max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                              ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                              number_iterations, iloop, feloop, &
                              converged, ires, cpu_finish - cpu_start, counter, &
                              minf,   &
                              x

                    call nlo_destroy(polish_opt)
                  enddo ! loop for two polishing algorithms
                else
                  polish_x = x
                  algorithm_polish = 'NONE'
                  call cpu_time(cpu_finish)
                  converged = 1
                  ! write results in file
                  write(3,*)  thread_no, nthreads, iseed, function_name, dim, &
                              algorithm, algorithm_polish, &
                              max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                              ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                              number_iterations, iloop, feloop, &
                              converged, ires, cpu_finish - cpu_start, counter, &
                              minf,   &
                              polish_x

                ! write interemediate result (before polishing part - only difference with above is last variable x vs polish_x)
                write(4,*)  thread_no, nthreads, iseed, function_name, dim, &
                            algorithm, algorithm_polish, &
                            max_evalstep, max_eval_mlsl_local, maxeval_polish, &
                            ftol_mlsl_local, ftol_polish, xtol_polish, x0(1), &
                            number_iterations, iloop, feloop, &
                            converged, ires, cpu_finish - cpu_start, counter, &
                            minf,   &
                            x


                endif ! end of local vs global algo loop

            endif ! if loop for convergence of optimization
		end do ! loop for different maxeval / feloop (run 50 times; here same random draw for each 50 times, before moving to next implementation)

    ! destroy the object opt
    call nlo_destroy(opt)
    if (algoloop == 3) then
      call nlo_destroy(local_opt)
		endif
    if (algoloop == 6) then
      call nlo_destroy(local_opt)
		endif
    if (algoloop == 7) then
      call nlo_destroy(local_opt)
		endif
    if (algoloop == 8) then
      call nlo_destroy(local_opt)
		endif

	enddo
  !$OMP END PARALLEL DO
enddo

print*, "--------------"
print*, "END OF PROGRAM"
contains

END PROGRAM main

!------------------------------------------------------------
!               INCOME PROCESS FUNCTION
!------------------------------------------------------------
subroutine incomeprocess(val, n, x, grad, need_gradient)

   USE OBJECTIVE, only : objFun !, dfovec, initial0 !,GetData
   use myparamsNLOPT, only: counter

IMPLICIT NONE

    double precision val, x(n), grad(n), xx(n)
    integer n, need_gradient, idim, igrad



    counter = counter + 1
    val = objFun(x)    ! from module objective
    write(*,*) val
end subroutine incomeprocess
