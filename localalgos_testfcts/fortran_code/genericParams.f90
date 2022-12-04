MODULE genericParams
    !====================================================================
    !Define variables we need for the global optimization.
    !Note, you should not put function specific parameters here. Use
    !a separate module.
    !====================================================================
    USE nrtype
    USE monteCarloParams
     
    IMPLICIT NONE

    !INTEGER(I4B) :: p_default_alg = 0   !default algorithm. 0-bobyqa, 1-amoeba, 2-dfpmin
    INTEGER(I4B) :: p_nx, p_nmom, p_iprint=2, p_maxeval, p_ninterppt ! pn_x is dimension of domain

    INTEGER(I4B):: p_qr_ndraw, p_maxpoints
    !INTEGER(I4B) :: p_searchType, p_lotteryPoints ! not used for local searches
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_wspace
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_init
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_range, p_bound
    !INTEGER(I4B), PARAMETER  :: p_nsmplx=100
        ! initialize the seed for random numbers
    !INTEGER(I4B), PARAMETER :: p_SEED=314159265
    !INTEGER(I4B), parameter :: sobdrop=0 !number of sobol points that are not kept (among first ones generated)
    !REAL(DP), parameter :: p_tolf_amoeba= 1.0d-3 ! Amoeba variables

    ! amoeba parameters
    REAL(DP), parameter :: p_tolf_amoeba= 1D-18 !0.000000000000000001D0
    INTEGER(I4B), parameter :: p_itmax_amoeba=1000 ! Amoeba variables

    ! dfpmin parameters
    REAL(DP), parameter :: dfpmin_tol = 1D-12 !0.000000000001D0

    ! dfnls parameters
    REAL(DP) :: itratio, rhobeg, rhoend

    !REAL(DP) :: global_min_f
    !REAL(DP), dimension(:), ALLOCATABLE :: global_min_x

! CONTAINS
!   ! needed for income_process
!   SUBROUTINE initialize() !(seqNu, update, config)
!       !Initialize the instance. This includes parsing the configuration file,
!       !setting parameter values, etc.
!
!       !INTEGER(I4B), INTENT(IN) :: seqNu
!       !LOGICAL, INTENT(IN), OPTIONAL :: update
!       !CHARACTER(LEN=25), INTENT(IN), OPTIONAL :: config
!       INTEGER(I4B) :: i, n, fileDesc, openStat
!       INTEGER, DIMENSION(:), ALLOCATABLE :: seed
!
!       call random_seed(size = p_nx) !n)
!       allocate(seed(p_nx)) !n))
!       seed(1)=123456
!       call random_seed(put = seed)
!       deallocate(seed)
!   END SUBROUTINE initialize

END MODULE genericParams
