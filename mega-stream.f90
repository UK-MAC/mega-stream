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
!
! This aims to test the theory that streaming many large arrays causes memory
! bandwidth limits not to be reached, and latency becomes a dominating factor.
! We run a kernel with a similar form to the original triad, but with more than
! 3 input arrays.
!
! The main kernel computes:
! r(i,j,k) = q(i,j,k) + a(i)*x(i,j) + b(i)*y(i,j) + c(i)*z(i,j)
! total(j,k) = SUM(r(:,j,k))
!

PROGRAM megastream

  USE omp_lib

  IMPLICIT NONE

!
! Arrays are defined in terms of 3 sizes
! The large arrays are of size SMALL*MEDIUM*LARGE and are indexed with 3 indicies.
! The medium arrays are of size SMALL*MEDIUM and are indexed with 2 indicies.
! The small arrays are of size SMALL and are indexed with 1 index.
!
!  By default the large array has 2^27 elements, and the small array has 64 elements (2^6).
!

  INTEGER :: L_size = 4096 ! 2^12
  INTEGER :: M_size =  512 ! 2^9
  INTEGER :: S_size =   64 ! 2^6
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: q, r
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: x, y, z
  REAL(8), DIMENSION(:), ALLOCATABLE :: a, b, c
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: total
  REAL(8) :: tmp

  INTEGER :: ntimes = 100
  REAL(8), DIMENSION(:), ALLOCATABLE :: timings
  REAL(8) :: tick, tock

  REAL(8) :: size

  REAL(8) :: gold, gold_total

  INTEGER :: i, j, k, t

! Tollerance with which to check final array values */
  REAL(8) :: TOLR = 1.0E-15

! Print information
  WRITE(*, *) 'MEGA-STREAM! - v0.2.1'
  WRITE(*, *)
  WRITE(*, '(a,I,a,f8.1,a)') 'Small arrays: ', S_size, ' elements', S_size*8*1.0E-3, 'KB'
  WRITE(*, '(a,I,a,I,a,f8.1,a)') 'Medium arrays: ', S_size, 'x', M_size, ' elements', S_size*M_size*8*1.0E-6, 'MB'
  WRITE(*, '(a,I,a,I,a,I,a,f8.1,a)') 'Large arrays: ', S_size, 'x', M_size, 'x', L_size, ' elements', S_size*M_size*L_size*8*1.0E-6, 'MB'
  WRITE(*, '(a,f8.1,a)') 'Memory footprint: ', 8 * 1.0E-6 * ( &
    2.0*L_size*M_size*S_size + &
    3.0*M_size*S_size +        &
    3.0*S_size +               &
    L_size*M_size)             &
    , ' MB'
  WRITE(*, *) 'Running ', ntimes, 'times'
  WRITE(*, *)

! Total memory moved - the arrays plus an extra total as update is +=
  size = 8 * 1.0E-6 * ( &
    2.0*L_size*M_size*S_size + &
    3.0*M_size*S_size +        &
    3.0*S_size +               &
    2.0*L_size*M_size)

! Allocate memory
  ALLOCATE(q(S_size, M_size, L_size))
  ALLOCATE(r(S_size, M_size, L_size))
  ALLOCATE(x(S_size, M_size))
  ALLOCATE(y(S_size, M_size))
  ALLOCATE(z(S_size, M_size))
  ALLOCATE(a(S_size))
  ALLOCATE(b(S_size))
  ALLOCATE(c(S_size))
  ALLOCATE(total(M_size, L_size))
  ALLOCATE(timings(ntimes))

! Initalise the data */
!$OMP PARALLEL

!$OMP DO
  DO k = 1, L_size
    DO j = 1, M_size
      DO i = 1, S_size
        q(i,j,k) = 0.1_8
        r(i,j,k) = 0.0_8
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO
  DO j = 1, M_size
    DO i = 1, S_size
      x(i,j) = 0.2_8
      y(i,j) = 0.3_8
      z(i,j) = 0.4_8
    END DO
  END DO
!$OMP END DO

!$OMP DO
  DO i = 1, S_size
    a(i) = 0.6_8
    b(i) = 0.7_8
    c(i) = 0.8_8
  END DO
!$OMP END DO

!$OMP DO
  DO k = 1, L_size
    DO j = 1, M_size
      total(j, k) = 0.0_8
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

! Run the kernel multiple times
  DO t = 1, ntimes
    tick = omp_get_wtime()

!**************************************************************************
!* Kernel
!*************************************************************************/
!$OMP PARALLEL DO PRIVATE(tmp)
    DO k = 1, L_size
     DO j = 1, M_size
        tmp = 0.0_8
        !$OMP SIMD REDUCTION(+:tmp)
        DO i = 1, S_size
          r(i,j,k) = q(i,j,k) + a(i)*x(i,j) + b(i)*y(i,j) + c(i)*z(i,j)
          tmp = tmp + r(i,j,k)
        END DO ! i
        total(j,k) = total(j,k) + tmp
      END DO ! j
    END DO ! k
!$OMP END PARALLEL DO

    tock = omp_get_wtime()
    timings(t) = tock-tick

  END DO ! t


! Check the results
  gold = 0.1 + 0.2*0.6 + 0.3*0.7 + 0.4*0.8
  gold_total = gold*S_size*ntimes
  IF (ANY(ABS(r - gold) > TOLR)) THEN
    WRITE(*, *) 'Results incorrect'
  END IF
  IF (ANY(ABS(total - gold_total) > TOLR)) THEN
    WRITE(*, *) 'Reduciton incorrect'
  END IF

! Print timings
  WRITE(*, '(a,a,a,a)') 'Bandwidth MB/s', 'Min time', 'Max time', 'Avg time'
  WRITE(*, '(f12.1,f11.6,f11.6,f11.6)') &
    size / MINVAL(timings(2:ntimes)), &
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

!
!void parse_args(int argc, char *argv[])
!{
!  for (int i = 1; i < argc; i++)
!  {
!    if (strcmp(argv[i], "--large") == 0)
!    {
!      L_size = atoi(argv[++i]);
!    }
!    else if (strcmp(argv[i], "--medium") == 0)
!    {
!      M_size = atoi(argv[++i]);
!    }
!    else if (strcmp(argv[i], "--small") == 0)
!    {
!      S_size = atoi(argv[++i]);
!    }
!    else if (strcmp(argv[i], "--swap") == 0)
!    {
!      int tmp = L_size;
!      L_size = M_size;
!      M_size = tmp;
!    }
!    else if (strcmp(argv[i], "--ntimes") == 0)
!    {
!      ntimes = atoi(argv[++i]);
!      if (ntimes < 2)
!      {
!        fprintf(stderr, "ntimes must be 2 or greater\n");
!        exit(EXIT_FAILURE);
!      }
!    }
!    else if (strcmp(argv[i], "--help") == 0)
!    {
!      printf("Usage: %s [OPTION]\n", argv[0]);
!      printf("\t --large n \tSet size of large dimension\n");
!      printf("\t --medium n \tSet size of medium dimension\n");
!      printf("\t --small n \tSet size of small dimension\n");
!      printf("\t --swap\tSwap medium and large sizes over\n");
!      printf("\t --ntimes n\tRun the benchmark n times\n");
!      printf("\n");
!      printf("\t Large  is %12d elements\n", LARGE);
!      printf("\t Medium is %12d elements\n", MEDIUM);
!      printf("\t Small  is %12d elements\n", SMALL);
!      exit(EXIT_SUCCESS);
!    }
!    else
!    {
!      fprintf(stderr, "Unrecognised argument \"%s\"\n", argv[i]);
!      exit(EXIT_FAILURE);
!    }
!  }
!}
