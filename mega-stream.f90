PROGRAM megastreamftn

  USE omp_lib

  IMPLICIT NONE
  
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: r, q
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: a, x

  REAL(8), ALLOCATABLE, DIMENSION(:) :: times

  REAL(8) :: size
 
  INTEGER :: L_size, S_size

  INTEGER :: i, t, ntimes

  DOUBLE PRECISION :: tick, tock

  CHARACTER(len=32) :: arg

  S_size = 128

! Check if --ssize is set

  IF (iargc() > 1) THEN
    CALL getarg(1, arg)
    IF (arg == "--ssize") THEN
      CALL getarg(2, arg)
      READ(arg,*) S_size
      PRINT *, "Setting S_size"
    ENDIF
  ENDIF

  L_size = 134217728 / S_size

  ntimes = 10

  size = 8.0 * (2.0*L_size*S_size + 2.0*S_size) * 1.0E-6
  
  ALLOCATE(r(S_size,L_size))
  ALLOCATE(q(S_size,L_size))
  ALLOCATE(a(S_size))
  ALLOCATE(x(S_size))
  ALLOCATE(times(ntimes))

!$OMP PARALLEL WORKSHARE
  r = 0.0_8
  q = 0.1_8
  x = 0.2_8
  a = 0.6_8
!$OMP END PARALLEL WORKSHARE

  iters: DO t = 1, ntimes

    tick = omp_get_wtime()

    !$OMP PARALLEL DO
      DO i = 1, L_size
        r(:,i) = q(:,i) + a(:)*x(:)
      END DO
    !$OMP END PARALLEL DO

    tock = omp_get_wtime()

    times(t) = tock-tick

  END DO iters

  print *, size/MINVAL(times(2:))
  
  DEALLOCATE(r)
  DEALLOCATE(q)
  DEALLOCATE(a)
  DEALLOCATE(x)
  DEALLOCATE(times)

END PROGRAM megastreamftn

