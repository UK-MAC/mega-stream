!
! Copyright 2016 Tom Deakin, University of Bristol
!
! This file is part of mega-stream.
!
! mega-stream is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! mega-stream is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with mega-stream.  If not, see <http://www.gnu.org/licenses/>.
!
! Test whether adding in a small array to a STREAM copy reduces memory bandwidth
!

PROGRAM megastreamftn

  USE omp_lib

  IMPLICIT NONE
  
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: r, q
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: a, b, c
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: x, y, z

  REAL(8), ALLOCATABLE, DIMENSION(:) :: times

  REAL(8) :: size
 
  INTEGER :: L_size, S_size

  INTEGER :: i, t, ntimes = 10

  DOUBLE PRECISION :: tick, tock

  CHARACTER(len=32) :: arg


! Set the default S_size
  S_size = 128

! --ssize option can be used to change S_size
  IF (iargc() > 1) THEN
    CALL getarg(1, arg)
    IF (arg == "--ssize") THEN
      CALL getarg(2, arg)
      READ(arg,*) S_size
      PRINT *, "Setting S_size"
    ENDIF
  ENDIF


! Set large size so total array is 2^27 doubles
  L_size = 134217728 / S_size


  ! Total memory movement in GB
  size = 8.0 * (2.0*L_size*S_size + 2.0*S_size) * 1.0E-6
  
  ALLOCATE(r(S_size,L_size))
  ALLOCATE(q(S_size,L_size))
  ALLOCATE(a(S_size))
  ALLOCATE(b(S_size))
  ALLOCATE(c(S_size))
  ALLOCATE(x(S_size))
  ALLOCATE(y(S_size))
  ALLOCATE(z(S_size))

  ALLOCATE(times(ntimes))


!Initilise data
!$OMP PARALLEL WORKSHARE
  r = 0.0_8
  q = 0.1_8
  x = 0.2_8
  y = 0.3_8
  z = 0.4_8
  a = 0.6_8
  b = 0.7_8
  c = 0.8_8
!$OMP END PARALLEL WORKSHARE

!-----------
! Main loop
!-----------

  iters: DO t = 1, ntimes

    tick = omp_get_wtime()

    !-----------
    ! Kernel
    !-----------
    !$OMP PARALLEL DO
    DO i = 1, L_size
      r(:,i) = q(:,i) + a(:)*x(:) + b(:)*y(:) + c(:)*z(:)
    END DO
    !$OMP END PARALLEL DO

    tock = omp_get_wtime()

    times(t) = tock-tick

  END DO iters


! Print bandwidth of fastest time, ignoring the first iteration
  print *, size/MINVAL(times(2:))
  
  DEALLOCATE(r)
  DEALLOCATE(q)
  DEALLOCATE(a)
  DEALLOCATE(b)
  DEALLOCATE(c)
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(times)

END PROGRAM megastreamftn

