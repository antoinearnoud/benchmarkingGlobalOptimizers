! -----------------------------------------------------------------------------------
!   Main Program:
!   This programm looks for the global minimum of test functions
!    using different algorithms  from NLOPT package and a Monte Carlo procedure.
!    The number of iterations (Monte Carlo) is defined by the number_iterations (we used 100).
!    The algorithms used are: CRS2_LM, ISRES, MLSL_LDS, STOGO, ESCH
!    Note: By default, STOGO is disactivated because it requires the installation
!    of the C++ extension of NLOPT package. To activate STOGO, set includeStogo = .TRUE.
!    The test functions used are levi13, griewank, rosenbrock, rastrigin.
!
!   This program uses NLOPT algorithms, courtesy of the NLOPT group at MIT.
!    see http://ab-initio.mit.edu/nlopt/
!
!   This program can be compiled with or without openmp for parallelization without
!    performance loss.
!    Openmp considerably decreases the time needed to run the program.
!
!   Program written for the paper Benchmarking Global Optimizers, 2020
!    by Antoine Arnoud, Fatih Guvenen and Tatjana Kleineberg
! ------------------------------------------------------------------------------------

program main
  !$ use omp_lib                                                                  ! omp_lib is used for parallelization with openmp
  use myparams_counter                                                            ! to get module variable counter
  use testFunctions
  implicit none
  ! functions variables
  !external                         :: levi13, griewank, rastrigin, rosenbrock	  ! test functions (if not using module testFunctions.f90)
  integer, dimension(2)            :: dim_list = (/ 2, 10 /)                     ! array of dimensions used
  integer                          :: dim                                        ! dimension of the problem
  double precision, allocatable		 :: x(:), polish_x(:), x0(:), x_random_mat(:,:) ! running point, polshing running point
  double precision, allocatable    :: solution_point(:)                          ! known solution of the minimization
  double precision                 :: fsol                                       ! known minimum
  double precision 				         :: minf                                       ! minimum value of objective function
  double precision, allocatable    :: lb(:), ub(:)                               ! lower and upper bounds (vectors) of the domain
  double precision 			           :: initpoint                                  !  initpoint: temporary varaible
  !NLopt algorithms variables
  integer*8                        :: opt, polish_opt, local_opt               ! pointers to NLopt algorithms
  character(35) 				           :: algorithm, algorithm_polish, function_name ! algorithm name, polishing phase algorithm name, function name
  integer                          :: thread_no, nthreads                      ! used with openmp
  real                             :: cpu_start, cpu_finish                    ! computes time needed for optimization
  integer, parameter               :: number_iterations = 100                  !number of implementations (Monte Carlo). used to compute performance profiles
  ! stopping criteria
  integer                          :: max_evalstep                             ! maximum number of function evaluations
  integer, parameter               :: max_eval_mlsl_local = 5000               ! maximum number of evaluation (stopping criterion) for local search.
  double precision                 :: ftol_mlsl_local                          ! f-f_{-1} < ftol_mlsl_local
  double precision, parameter      :: ftol_mlsl_local_3 = 1.D-3
  double precision, parameter      :: ftol_mlsl_local_8 = 1.D-8
  ! stopping criteria (polishing phase)
  integer                          :: maxeval_polish                           ! for example max_evalstep / 5 or max_evalstep/2
  double precision, parameter      :: ftol_polish = 1.D-8
  double precision, parameter      :: xtol_polish = 1.D-8
  integer 					               :: ires                                     ! success value
  ! other parameters
  integer, parameter               :: iseed = 123456
  logical, parameter               :: includeStogo = .TRUE. ! set to FALSE if don't want to run STOGO
  ! loops variables
  integer 					               :: iloop, algoloop, funcloop, dim_loop, feloop, polishloop, idim, jj ! loop varaibles
  integer            					     :: converged

  include 'nlopt.f'                               ! this is needed to call NLOPT optimizers

  ! create the file for output (add a space at the end of each variable name)
  open(3, file='nlopt_testfcts.txt')
  write(3,*)    'thread_no ', 'number_of_threads ', 'seed ', 'test_function ', 'dimension ', &
                'algorithm ', 'algorithm_polish ',  &
                'maxeval ', 'maxeval_mlsl_local ',  'maxeval_polish ', &
                'ftol_mlsl_local ', 'ftol_polish ', 'xtol_polish ' , 'startpt ', &
                'number_iterations ', 'iloop ', 'feloop ', &
                'converged ', 'ires ','cpu_time ', 'counter ',  &
                'minf ',  &
                'x1 ', 'x2 ', 'x3 ', 'x4 ', 'x5 ', 'x6 ', 'x7 ', 'x8 ', 'x9 ', 'x10 '

  ! saving results beefore polishing phase (in case need to re-run only the polishing phase)
  open(4, file='nlopt_testfcts_pre_polish.txt')
  write(4,*)    'thread_no ', 'number_of_threads ', 'seed ', 'test_function ', 'dimension ', &
                'algorithm ', 'algorithm_polish ',  &
                'maxeval ', 'maxeval_mlsl_local ',  'maxeval_polish ', &
                'ftol_mlsl_local ', 'ftol_polish ', 'xtol_polish ' , 'startpt ', &
                'number_iterations ', 'iloop ', 'feloop ', &
                'converged ', 'ires ', 'cpu_time ', 'counter ',  &
                'minf ',  &
                'x_prepolish_1 ', 'x_prepolish_2 ', 'x_prepolish_3 ', 'x_prepolish_4 ', 'x_prepolish_5 ', &
                'x_prepolish_6 ', 'x_prepolish_7 ', 'x_prepolish_8 ', 'x_prepolish_9 ', 'x_prepolish_10 '

  ! loop over dimensionality
  do dim_loop = 1, 2
    dim = dim_list(dim_loop)

    ! allocate matrices
    allocate(lb(dim), ub(dim), solution_point(dim))
    allocate(x(dim), polish_x(dim), x0(dim), x_random_mat(dim, number_iterations))

    ! loop over test functions
    do funcloop = 1, 4
      !read in the starting points
      if (funcloop==1) then
        function_name = "levi13"
        open(222, file= "init_points/init_points_levi.dat")
        fsol = 1.0                                ! minimum value if 0.0 but we added 1 to the function value
        do idim=1, dim
          solution_point(idim) = 1.0              ! known solution (used to evaluate success)
          lb(idim) = -10.0                        ! lower and upper bounds of the search domain
          ub(idim) = 10.0
        enddo
      else if (funcloop==2) then
        function_name = "griewank"
        open(222, file= "init_points/init_points_griewank.dat")
        ! do iloop = 1, number_iterations
        !   do jj = 1, dim
        !    read(221,*) initpoint
        !    x_random_mat(jj,iloop) = initpoint
        !   enddo
        ! enddo
        ! close(221)
        fsol = 1.0                                  ! minimum value if 0.0 but we added 1 to the function value
        do idim=1, dim
          solution_point(idim) = 0.0
          lb(idim) = -100.0
          ub(idim) = 95.0                         ! upper bound is 100 but we rescale the bounds asymetrically so that algorithms do not directly guess the right point by change (center of the search space)
        enddo
      else if (funcloop==3) then
        function_name = "rosenbrock"
        open(222, file= "init_points/init_points_rosen.dat")
        ! do iloop = 1, number_iterations
        !   do jj = 1, dim
        !    read(223,*) initpoint
        !    x_random_mat(jj,iloop) = initpoint
        !   enddo
        ! enddo
        ! close(223)
        fsol = 1.0
        do idim=1, dim
          solution_point(idim)=1.0
          lb(idim) = -100.0
          ub(idim) = 100.0
        enddo
      else
        function_name = "rastrigin"
        open(222, file= "init_points/init_points_rastrigin.dat")
        ! do iloop = 1, number_iterations
        !   do jj = 1, dim
        !    read(224,*) initpoint
        !    x_random_mat(jj,iloop) = initpoint
        !   enddo
        ! enddo
        ! close(224)
        fsol = 1.0
        do idim=1, dim
          solution_point(idim)=0.0
          lb(idim) = -5.12
          ub(idim) = 5.0                          ! upper bound is 5.12 but rescale the bounds asymetrically  so that id does not directly guess the right point
        enddo
      endif

      do iloop = 1, number_iterations
        do jj = 1, dim
         read(222,*) initpoint
         x_random_mat(jj,iloop) = initpoint
        enddo
      enddo
      close(222)

      ! loop over the different search algorithms
      do algoloop = 1, 11
        if ((includeStogo .eqv. .FALSE.) .and. (algoloop == 4)) then ! skip STOGO if includeStogo parameter is set to FALSE
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
          algorithm = 'NLOPT_G_MLSL_LDS_NELDERMEAD_8'
        else if (algoloop == 7) then
          algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_3'
        else if (algoloop == 8) then
          algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_8'
        else if (algoloop == 9) then
          algorithm = 'NLOPT_LD_LBFGS'              !similar to dfpmin
        else if (algoloop == 10) then
          algorithm = 'NLOPT_LN_NELDERMEAD'
        else if (algoloop == 11) then
          algorithm = 'NLOPT_LN_BOBYQA'
        endif

        ! loop over iterations (each starts with different initial point)
        ! note: number_iterations should not be in a private() or shared() statement since it is a parameter
        ! old one, delete: ! $omp parallel do if(number_iterations.gt.10 .AND. algoloop /= 3) default(none), &

        !OMP PARALELLIZATION STARTS HERE
        !declare paralellization parameters. uses openmp only if more than 10 iterations.
        !$omp parallel do if(number_iterations.gt.10) default(none), &
        !$omp& private(cpu_start, cpu_finish), &
        !$OMP& private(x0, x, polish_x, minf, ires, feloop, polishloop),  &
        !$OMP& private(max_evalstep, maxeval_polish, ftol_mlsl_local), &
        !$OMP& private(opt, polish_opt, local_opt, algorithm_polish, nthreads, thread_no, converged), &
        !$OMP& shared(x_random_mat, algoloop, dim, lb, ub, funcloop), &
        !$OMP& shared(function_name, algorithm, fsol, solution_point)
        do iloop = 1, number_iterations
          !get number of threads available
          nthreads = -1 ! USED WHEN NOT USING OMP; overwritten by line below when using OMP
          !$ nthreads = OMP_GET_NUM_THREADS()
          thread_no = -1 ! used when not using OMP; overwritten by lines below when using OMP
          !$ thread_no = OMP_GET_THREAD_NUM()
          ! $ print*, "number of threads: ", nthreads
          ! $ print*, "thread numero: " , thread_no

          ! Note: delete opt object at the end of this loop (and re-create)
          !    because dim (dimension) and algo (algorithm) in opt objects are immmutable.

          ! set the seed for each run
          call nlosr(iseed)

          ! value written when not using it (MLSL)
          ftol_mlsl_local = 0.0

          ! creates the opt object (indicates the algorithm used and the dimensionality)
          if (algoloop == 1) then
            !algorithm = 'NLOPT_GN_CRS2_LM'
            call nlo_create(opt, NLOPT_GN_CRS2_LM, dim)
          else if (algoloop == 2) then
            !algorithm = 'NLOPT_GN_ISRES'
            call nlo_create(opt, NLOPT_GN_ISRES, dim)
          else if (algoloop == 3) then
            !algorithm = 'NLOPT_G_MLSL_LDS_3'
            call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
            ! create the local opt (MLSL needs a local optimizer. we use Nelder Mead here, and BOBYQA below)
            call nlo_create(local_opt, NLOPT_LN_NELDERMEAD, dim)
            ftol_mlsl_local = ftol_mlsl_local_3
            call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_3)  ! it is ok to start with a high tolerance for local search and do a final local search at the end.
                                                              ! MLSL defaults to ftol_rel=10−15 and xtol_rel=10−7 for the local searches
            call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
            call nlo_set_local_optimizer(ires, opt, local_opt)
          else if (algoloop == 4) then
            !algorithm = 'NLOPT_GD_STOGO'
            call nlo_create(opt, NLOPT_GD_STOGO, dim)
          else if (algoloop == 5) then
            !algorithm = 'NLOPT_GN_ESCH'
            call nlo_create(opt, NLOPT_GN_ESCH, dim)
          else if (algoloop == 6) then
            !algorithm = 'NLOPT_G_MLSL_LDS_8'
            call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
            ! create the local opt (MLSL needs a local optimizer. we use Nelder Mead)
            call nlo_create(local_opt, NLOPT_LN_NELDERMEAD, dim)
            ftol_mlsl_local = ftol_mlsl_local_8
            call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_8)
            call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
            call nlo_set_local_optimizer(ires, opt, local_opt)
          else if (algoloop == 7) then ! MLSLS with BOBYQA local searches
            !algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_3'
            call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
            ! create the local opt (MLSL needs a local optimizer. we use BOBYQA)
            call nlo_create(local_opt, NLOPT_LN_BOBYQA, dim)
            ftol_mlsl_local = ftol_mlsl_local_3
            call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_3)
            call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
            call nlo_set_local_optimizer(ires, opt, local_opt)
          else if (algoloop == 8) then
            !algorithm = 'NLOPT_G_MLSL_LDS_BOBYQA_8'
            call nlo_create(opt, NLOPT_G_MLSL_LDS, dim)
            call nlo_create(local_opt, NLOPT_LN_BOBYQA, dim)
            ftol_mlsl_local = ftol_mlsl_local_8
            call nlo_set_ftol_rel(ires, local_opt, ftol_mlsl_local_8)
            call nlo_set_maxeval(ires, local_opt, max_eval_mlsl_local)
            call nlo_set_local_optimizer(ires, opt, local_opt)
          else if (algoloop == 9 ) then
            !algorithm = 'NLOPT_LD_LBFGS'
            call nlo_create(opt, NLOPT_LD_LBFGS, dim)
          else if (algoloop == 10) then
            !algorithm = 'NLOPT_LN_NELDERMEAD'
            call nlo_create(opt, NLOPT_LN_NELDERMEAD, dim)
          else if (algoloop == 11) then
            !algorithm = 'NLOPT_LN_BOBYQA'
            call nlo_create(opt, NLOPT_LN_BOBYQA, dim)
          end if

          ! set lower and upper bounds
          call nlo_set_lower_bounds(ires, opt, lb)
          call nlo_set_upper_bounds(ires, opt, ub)
          !other stopping criteria to avoid early-stops of algos.
          !these are hurting the chances of NLOPT.
          !call nlo_set_ftol_rel(ires, opt, 1.D-20)
          !call nlo_set_xtol_rel(ires, opt, 1.D-20)
          !call nlo_set_ftol_abs(ires, opt, 1.D-20)
          !call nlo_set_xtol_abs(ires, opt, 1.D-20)

          ! set the objective function to minimize
          ! last argument = 0 means no additional data to be passed to the function
          if (funcloop == 1) then
            call nlo_set_min_objective(ires, opt, levi13, 0)
          else if (funcloop == 2) then
            call nlo_set_min_objective(ires, opt, griewank, 0)
          else if (funcloop == 3) then
            call nlo_set_min_objective(ires, opt, rosenbrock, 0)
          else if (funcloop == 4) then
            call nlo_set_min_objective(ires, opt, rastrigin, 0)
  	      else
    		    print*,"error: unknown function"
            stop                                    ! stops the program if unknown function
          endif

          call cpu_time(cpu_start)

          ! random initial guess point (parallelization along iloop)
          x0 = x_random_mat(:,iloop)
          max_evalstep = 0          ! initiliaze. used in the first loop below.
          do feloop = 1, 30         ! FIXME: use an array instead of looping?
            x = x0       						!reset x to starting point (before this, x is the output from last optimization)
            if (feloop <= 9) then
              max_evalstep = max_evalstep + 50*dim               !feloop * 50*dim
            endif
            if (feloop >= 10) then
              max_evalstep = max_evalstep + dim*10000/20          !feloop *dim*10000/20   !1000=increments for 2 dim; 5000 increments for 10 dim
            endif

            call nlo_set_maxeval(ires, opt, max_evalstep)

            ! initialize counter of number of calls to the function during search
            counter = 0

            ! find optimum
            call nlo_optimize(ires, opt, x, minf) 		! optimizes the object opt. returns x and minf. ires is message or error.

            ! used if no convergence. if convergence, overwritten below.
            call cpu_time(cpu_finish)

            !print results of search (if ires<0: issue, -4: Halted because roundoff errors limited progress)
            if (ires.lt.0 .AND. ires /= -4) then
              ! write(*,*) 'nlopt failed! error code: ', ires
              ! write(*,*) ' algorithm ', algorithm
              ! write(*,*) ' function ', function_name
              ! write(*,*) ' dimension ' , dim
              ! write(*,*) ' maxeval ', max_evalstep
              ! write(*,*) 'implementation number: ', iloop
              ! write(*,*) 'initial point', x
              converged = -1
              algorithm_polish = 'NONE'
              maxeval_polish = -999999
              minf = -999999
              if (dim == 2) then
                polish_x = (/ -999999, -999999 /)
              endif
              if (dim == 10) then
                polish_x = (/ -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999 /)
              endif
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
            else
              !write(*,*) 'algorithm: ', algorithm
              !write(*,*) 'function: ' , function_name
              !write(*,*) 'dimension: ', dim
              !write(*,*) 'maxevalstep: ', max_evalstep
              !write(*,*) 'iteration: ', iloop
              !write(*,*) 'x0 = ', x0
              !write(*,*) 'found min at ', x
              !write(*,*) 'min val = ', minf
              !write(*,*) 'found after ', counter, 'iterations'
              !write(*,*) 'found in ', cpu_finish - cpu_start, 'cpu seconds'
              !write(*,*) 'ires code: ', ires
              !write(*,*) 'real minval at: ', solution_point
              !write(*,*) 'real minval = ', fsol
              !write(*,*) '------------------------------------------------------------------------------'
              !write(*,*) ''

              ! set the polishing optimization
              maxeval_polish = max_evalstep/5 ! integer / integer -> integer
              !Call the local optimizer for polishing part
              !not for NLOPT_LD_LBFGS (=9), NELDER-MEAD (=10), BOBYQA (=11)
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

                  if (funcloop == 1) then
                    call nlo_set_min_objective(ires, polish_opt, levi13, 0)
                  else if (funcloop == 2) then
                    call nlo_set_min_objective(ires, polish_opt, griewank, 0)
                  else if (funcloop == 3) then
                    call nlo_set_min_objective(ires, polish_opt, rosenbrock, 0)
                  else if (funcloop == 4) then
                    call nlo_set_min_objective(ires, polish_opt, rastrigin, 0)
                  else
                    print*, "error: unknown function"
                    stop
                  endif

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
            endif ! if-loop for convergence of the global optmizer
          end do ! loop for different maxeval / feloop (run 30 steps of performance profile with same random draw, before moving to next implementation)
          ! need to destroy the opt *before* getting out of the parallel loop
          ! because they are private so once out of OMP loop they are not initialized / defined anymore
          ! so destroy(opt) would return an error
          call nlo_destroy(opt)
          if ((algoloop == 3) .OR. (algoloop == 6) .OR. (algoloop == 7) .OR. (algoloop == 8)) then
            call nlo_destroy(local_opt)
          endif
        end do	! end loop over implementations (e.g. new random draw)
        !$OMP END PARALLEL DO ! used for parallalization with openmp
      enddo 		! endo algoloop
    enddo 		! end loop on functions
    deallocate(x, polish_x, x0, x_random_mat, lb, ub, solution_point)
  enddo 			! end loop on dimension
  close(3)    ! close file nlopt_testfcts.dat
  close(4)    ! close file nlopt_testfcts.dat
End program main

!------------------------------------------------------------
!               TEST FUNCTIONS
!include 'testFunctions.f90'
!------------------------------------------------------------

! END OF FILE.
