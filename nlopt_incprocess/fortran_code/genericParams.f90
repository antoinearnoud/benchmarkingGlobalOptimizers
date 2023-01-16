!=========================================================================
! 2. GENERICPARAMS: Define variables we need for the global optimization.
    !Note, you should not put function specific parameters here. Use
    !a separate module.
!=========================================================================
MODULE genericParams
USE nrtype
!USE stateControl
IMPLICIT NONE

    INTEGER(I4B) :: p_default_alg = 0 !default algorithm. 0-bobyqa, 1-amoeba, 2-dfpmin
    INTEGER(I4B) :: p_nx, p_nmom, p_iprint=2, p_maxeval, p_ninterppt
    INTEGER(I4B):: p_qr_ndraw, p_maxpoints
    INTEGER(I4B) :: p_searchType, p_lotteryPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_wspace
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_init
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_range, p_bound
    INTEGER(I4B), PARAMETER  :: p_nsmplx=100
    ! initialize the seed for random numbers
    INTEGER(I4B), PARAMETER :: p_SEED=314159265
    REAL(DP), parameter :: p_tolf_amoeba= 1.0d-8 ! Amoeba variables
    INTEGER(I4B), parameter :: p_itmax_amoeba=1000 ! Amoeba variables

    INTEGER(I4B) :: fe_counter ! = 0 no need to initialize
    REAL(DP) :: cpuTime = 0

    ! declaration used for openmp (need to declare threadprivate variables in their module)
    ! !$omp threadprivate()


CONTAINS
  SUBROUTINE initialize(seqNu, update, config)
      !Initialize the instance. This includes parsing the configuration file,
      !setting parameter values, etc.

      INTEGER(I4B), INTENT(IN) :: seqNu
      LOGICAL, INTENT(IN), OPTIONAL :: update
      CHARACTER(LEN=25), INTENT(IN), OPTIONAL :: config
      INTEGER(I4B) :: i, n, fileDesc, openStat
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      call random_seed(size = n)
      allocate(seed(n))
      seed(1)=123456
      call random_seed(put = seed)
      deallocate(seed)
  END SUBROUTINE initialize


END MODULE genericParams
