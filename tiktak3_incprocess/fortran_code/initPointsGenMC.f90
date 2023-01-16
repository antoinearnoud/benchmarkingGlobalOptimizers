PROGRAM initPointsGenMC
  ! This program uses config.txt to generates 1000 random points
  ! in the domain space defined in config.txt
  ! The points are draw uniformely on the p_nx dimensional domain space.
  ! The points are used for Monte Carlo simulation.

  IMPLICIT NONE

  CHARACTER(LEN=100) 	:: line
  REAL			:: mynumber
  INTEGER			:: p_nx, nmom, maxeval, qr_ndraw, maxpoints, searchType, lotteryPoints
  REAL, ALLOCATABLE	:: p_bound(:,:)
  LOGICAL			:: parsedLine
  INTEGER			:: i, jj, k


  ! Read dimension and bounds from configfile
  OPEN(UNIT=222, FILE="config.txt", STATUS='old', ACTION='read')

  !Parse the basic parameters in the first line
  parsedLine = .FALSE.
  DO WHILE (parsedLine .EQV. .FALSE.)
     READ(222,'(A100)') line

     IF (line(1:1)=='!') THEN
        CYCLE
     END IF

     parsedLine = .TRUE.
     READ(line,*) p_nx, nmom, maxeval, qr_ndraw, maxpoints, searchType, lotteryPoints
  END DO

  ALLOCATE(p_bound(p_nx,2))

  !Parse the bounds
  parsedLine = .FALSE.
  DO WHILE (parsedLine .EQV. .FALSE.)
     READ(222,'(A100)') line

     IF (line(1:1)=='!') THEN
        CYCLE
     END IF

     parsedLine = .TRUE.

     READ(line,*) p_bound(1,1), p_bound(1,2)

     DO i=2,p_nx
        READ(222,*) p_bound(i,1), p_bound(i,2)
     END DO
  END DO

  ! Generate 1000 initial points scaled to the bounds
  DO jj = 1,1000
     DO i = 1,p_nx
        ! generate a random number between 0 and 1
        CALL RANDOM_NUMBER(mynumber)        	
        ! scale the number to the bounds
        mynumber = mynumber * (p_bound(i,2) - p_bound(i,1)) + p_bound(i,1)
        ! write the scaled number in dat file
        OPEN(unit = 111, file = "init_points.dat", position = "append", STATUS='unknown')
        WRITE(111,*) mynumber
        CLOSE(111)
        CLOSE(222)
     ENDDO
  ENDDO
END PROGRAM initPointsGenMC
