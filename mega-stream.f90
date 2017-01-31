!
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
! This aims to investigate the limiting factor for a simple kernel, in particular
! where bandwidth limits not to be reached, and latency becomes a dominating factor.
!

PROGRAM megastream

  USE omp_lib

  IMPLICIT NONE

  ! Constant parameters

  ! Arrays are defined in terms of 3 sizes: inner, middle and outer.
  ! The large arrays are of size inner*middle*middle*middle*outer.
  ! The medium arrays are of size inner*middle*middle*outer.
  ! The small arrays are of size inner and are indexed with 1 index.

  INTEGER, PARAMETER :: OUTER = 64
  INTEGER, PARAMETER :: MIDDLE = 16
  INTEGER, PARAMETER :: INNER = 128

  ! Tollerance with which to check final array values
  REAL(8), PARAMETER :: TOLR = 1.0E-15

  ! Starting values
  REAL(8), PARAMETER :: R_START = 0.0_8
  REAL(8), PARAMETER :: Q_START = 0.01_8
  REAL(8), PARAMETER :: X_START = 0.02_8
  REAL(8), PARAMETER :: Y_START = 0.03_8
  REAL(8), PARAMETER :: Z_START = 0.04_8
  REAL(8), PARAMETER :: A_START = 0.06_8
  REAL(8), PARAMETER :: B_START = 0.07_8
  REAL(8), PARAMETER :: C_START = 0.08_8

  ! Variables

  ! Default strides
  INTEGER :: Ni = INNER
  INTEGER :: Nj = MIDDLE
  INTEGER :: Nk = MIDDLE
  INTEGER :: Nl = MIDDLE
  INTEGER :: Nm = OUTER

  ! Number of iterations to run benchmark
  INTEGER :: ntimes = 100

  ! Arrays
  REAL(8), DIMENSION(:,:,:,:,:), POINTER :: q, r, ptr_tmp
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: x, y, z
  REAL(8), DIMENSION(:), ALLOCATABLE :: a, b, c
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: total

  REAL(8), DIMENSION(:), ALLOCATABLE :: timings
  REAL(8) :: tick, tock

  REAL(8) :: moved

  INTEGER :: i, j, k, l, m, t

  ! Print information
  WRITE(*, *) 'MEGA-STREAM! - v0.3'
  WRITE(*, *)
  !WRITE(*, '(a,I,a,f8.1,a)') 'Small arrays: ', S_size, ' elements', S_size*8*1.0E-3, 'KB'
  !WRITE(*, '(a,I,a,I,a,f8.1,a)') 'Medium arrays: ', S_size, 'x', M_size, ' elements', S_size*M_size*8*1.0E-6, 'MB'
  !WRITE(*, '(a,I,a,I,a,I,a,f8.1,a)') 'Large arrays: ', S_size, 'x', M_size, 'x', L_size, ' elements', S_size*M_size*L_size*8*1.0E-6, 'MB'
  !WRITE(*, '(a,f8.1,a)') 'Memory footprint: ', 8 * 1.0E-6 * ( &
    !2.0*L_size*M_size*S_size + &
    !3.0*M_size*S_size +        &
    !3.0*S_size +               &
    !L_size*M_size)             &
    !, ' MB'
  WRITE(*, *) 'Running ', ntimes, 'times'
  WRITE(*, *)


  ! Allocate memory
  NULLIFY(q, r)
  ALLOCATE(q(Ni,Nj,Nk,Nl,Nm))
  ALLOCATE(r(Ni,Nj,Nk,Nl,Nm))
  ALLOCATE(x(Ni,Nj,Nk,Nm))
  ALLOCATE(y(Ni,Nj,Nl,Nm))
  ALLOCATE(z(Ni,Nk,Nl,Nm))
  ALLOCATE(a(Ni))
  ALLOCATE(b(Ni))
  ALLOCATE(c(Ni))
  ALLOCATE(total(Nj,Nk,Nl,Nm))

  ALLOCATE(timings(ntimes))

  ! Initalise the data */
!$OMP PARALLEL
  ! q and r
!$OMP DO
  DO m = 1, Nm
    DO l = 1, Nl
      DO k = 1, Nk
        DO j = 1, Nj
          DO i = 1, Ni
            r(i,j,k,l,m) = R_START
            q(i,j,k,l,m) = Q_START
          END DO
        END DO
       END DO
     END DO
   END DO
!$OMP END DO

  ! x
!$OMP DO
  DO m = 1, Nm
    DO k = 1, Nk
      DO j = 1, Nj
        DO i = 1, Ni
          x(i,j,k,m) = X_START
        END DO
      END DO
     END DO
   END DO
!$OMP END DO

  ! y
!$OMP DO
  DO m = 1, Nm
    DO l = 1, Nl
      DO j = 1, Nj
        DO i = 1, Ni
          y(i,j,l,m) = Y_START
        END DO
      END DO
     END DO
   END DO
!$OMP END DO

  ! z
!$OMP DO
  DO m = 1, Nm
    DO l = 1, Nl
      DO k = 1, Nk
        DO i = 1, Ni
          z(i,k,l,m) = Z_START
        END DO
      END DO
     END DO
   END DO
!$OMP END DO

  ! a, b, and c
!$OMP DO
  DO i = 1, Ni
    a(i) = A_START
    b(i) = B_START
    c(i) = C_START
  END DO
!$OMP END DO

  ! sum
!$OMP DO
  DO m = 1, Nm
    DO l = 1, Nl
      DO k = 1, Nk
        DO j = 1, Nj
          total(j,k,l,m) = 1.0_8
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

  ! Run the kernel multiple times
  DO t = 1, ntimes
    tick = omp_get_wtime()

    CALL kernel(Ni, Nj, Nk, Nl, Nm, r, q, x, y, z, a, b, c, total)

    ! Swap the pointers
    ptr_tmp => q
    q => r
    r => ptr_tmp
    NULLIFY(ptr_tmp)

    tock = omp_get_wtime()
    timings(t) = tock-tick

  END DO ! t


  ! Check the results
  WRITE(*, *) "Sum total: ", SUM(total)

  ! Print timings
  WRITE(*, '(a,a,a,a)') 'Bandwidth MB/s', 'Min time', 'Max time', 'Avg time'
  WRITE(*, '(f12.1,f11.6,f11.6,f11.6)') &
    moved / MINVAL(timings(2:ntimes)), &
    MINVAL(timings(2:ntimes)), &
    MAXVAL(timings(2:ntimes)), &
    SUM(timings(2:ntimes)) / (ntimes - 1)

! Deallocate memory
  DEALLOCATE(q, r)
  DEALLOCATE(a, b, c)
  DEALLOCATE(x, y, z)
  DEALLOCATE(total)
  DEALLOCATE(timings)

END PROGRAM megastream

!**************************************************************************
!* Kernel
!*************************************************************************/
SUBROUTINE kernel(Ni, Nj, Nk, Nl, Nm, r, q, x, y, z, a, b, c, total)

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: Ni, Nj, Nk, Nl, Nm

  REAL(8), DIMENSION(Ni, Nj, Nk, Nl, Nm), INTENT(INOUT) :: q, r
  REAL(8), DIMENSION(Ni, Nj, Nk, Nm), INTENT(INOUT) :: x
  REAL(8), DIMENSION(Ni, Nj, Nl, Nm), INTENT(INOUT) :: y
  REAL(8), DIMENSION(Ni, Nk, Nl, Nm), INTENT(INOUT) :: z
  REAL(8), DIMENSION(Ni), INTENT(IN) :: a, b, c
  REAL(8), DIMENSION(Nj, Nk, Nl, Nm), INTENT(INOUT) :: total

  ! Local variables
  INTEGER :: i, j, k, l, m
  REAL(8) :: tmp_r, tmp_total

!$OMP PARALLEL DO PRIVATE(tmp_r, tmp_total)
  DO m = 1, Nm
    DO l = 1, Nl
      DO k = 1, Nk
        DO j = 1, Nj
          tmp_total = 0.0_8
          !$OMP SIMD REDUCTION(+:tmp_total)
          DO i = 1, Ni
            ! Set r
            tmp_r = q(i,j,k,l,m) +  &
              a(i)*x(i,j,k,m) +     &
              b(i)*y(i,j,l,m) +     &
              c(i)*z(i,k,l,m)

            ! Update x, y and z
            x(i,j,k,m) = 0.2_8*tmp_r - x(i,j,k,m)
            y(i,j,l,m) = 0.2_8*tmp_r - y(i,j,l,m)
            z(i,k,l,m) = 0.2_8*tmp_r - z(i,k,l,m)

            ! Reduce over Ni
            tmp_total = tmp_total + tmp_r

            ! Save r
            r(i,j,k,l,m) = tmp_r

          END DO ! i

          total(j,k,l,m) = total(j,k,l,m) + tmp_total

        END DO ! j
      END DO ! k
    END DO ! l
  END DO ! m
!$OMP END PARALLEL DO

END SUBROUTINE kernel

