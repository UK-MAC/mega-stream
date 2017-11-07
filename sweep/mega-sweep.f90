!
!
! Copyright 2017 Tom Deakin, University of Bristol
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
program megasweep

  ! MEGA-SWEEP adds a KBA sweep along with the MPI communication to the MEGA-STREAM kernel

  use comms
  use omp_lib

  implicit none

  ! MPI variables
  integer :: rank, nprocs
  integer :: lrank, rrank ! neighbour ranks
  real(kind=8), dimension(:,:,:), allocatable :: buf ! comms buffer

  ! OpenMP variables
  integer :: nthreads

  ! Problem sizes
  integer :: nang, ng ! angles and groups
  integer :: nx, ny   ! global mesh size
  integer :: lnx, lny ! local mesh size
  integer :: chunk    ! chunk size
  integer :: nsweeps  ! sweep direction
  integer :: ntimes   ! number of times

  ! Arrays
  real(kind=8), dimension(:,:,:,:,:), pointer :: aflux0, aflux1 ! angular flux
  real(kind=8), dimension(:,:,:,:,:), pointer :: aflux_ptr      ! for pointer swap
  real(kind=8), dimension(:,:,:), allocatable :: sflux          ! scalar flux
  real(kind=8), dimension(:), allocatable :: mu, eta            ! angular cosines
  real(kind=8), dimension(:,:,:), allocatable :: psii, psij     ! edge angular flux
  real(kind=8), dimension(:), allocatable :: w                  ! scalar flux weights
  real(kind=8) :: v                                             ! time constant
  real(kind=8) :: dx, dy                                        ! cell size
  real(kind=8), dimension(:), allocatable :: pop                ! scalar flux reduction
  real(kind=8) :: total_pop

  ! Timers
  real(kind=8) :: start_time, end_time
  real(kind=8) :: total_time
  real(kind=8) :: timer
  real(kind=8), dimension(:), allocatable :: time

  ! Local variables
  integer :: t, g
  logical :: ydecomp = .false.
  real(kind=8) :: moved ! model of data movement
  real(kind=8) :: lmoved ! model of data movement for single MPI rank
  real(kind=8) :: tmoved ! model of data movement for single thread of single MPI rank

  call comms_init

  call comms_rank(rank)
  call comms_size(nprocs)

  !$omp parallel
    nthreads = omp_get_num_threads()
  !$omp end parallel

  ! Set default problem sizes
  nx = 128
  ny = 128
  ng = 16
  nang = 16
  nsweeps = 4
  chunk = 1
  ntimes = 500

  ! Read in command line arguments
  call parse_args(rank,nang,nx,ny,ng,chunk,ntimes,ydecomp)

  ! Check ny is split into even number of chunks
  if (.not. ydecomp .and. mod(ny,chunk) .ne. 0) then
    if (rank .eq. 0) then
      print *, "Chunk size must divide ny"
    end if
    stop
  end if

  if (ydecomp) then
    ! Decompose in y-dimension
    if (nprocs .gt. ny) then
      if (rank .eq. 0) then
        print *, "Too many processors to decompose ny"
      end if
      stop
    end if

    lnx = nx
    lny = ny / nprocs

    ! Share remainder cells for uneven decomposition
    ! Allows for flexible process counts
    if (mod(ny,nprocs) .ne. 0) then
      if (rank .lt. mod(ny,nprocs)) then
        lny = lny + 1
      end if
    end if
  else
    ! Decompose in x-dimension
    if (nprocs .gt. nx) then
      if (rank .eq. 0) then
        print *, "Too many processors to decompose nx"
      end if
      stop
    end if

    lnx = nx / nprocs
    lny = ny

    ! Share remainder cells for uneven decomposition
    ! Allows for flexible process counts
    if (mod(nx,nprocs) .ne. 0) then
      if (rank .lt. mod(nx,nprocs)) then
        lnx = lnx + 1
      end if
    end if
  end if
    
  ! Set neighbour ranks
  if (rank .ne. 0) then
    lrank = rank - 1
  else
    lrank = MPI_PROC_NULL
  end if
  if (rank .ne. nprocs-1) then
    rrank = rank + 1
  else
    rrank = MPI_PROC_NULL
  end if

  ! Allocate data
  allocate(aflux0(nang,lnx,lny,nsweeps,ng))
  allocate(aflux1(nang,lnx,lny,nsweeps,ng))
  allocate(sflux(lnx,lny,ng))
  allocate(mu(nang))
  allocate(eta(nang))
  if (ydecomp) then
    allocate(psii(nang,lny,ng))
    allocate(psij(nang,chunk,ng))
  else
    allocate(psii(nang,chunk,ng))
    allocate(psij(nang,lnx,ng))
  end if
  allocate(w(nang))
  allocate(pop(ng))
  allocate(buf(nang,chunk,ng))

  ! Initilise data
  !$omp parallel do
  do g = 1, ng
    aflux0(:,:,:,:,g) = 1.0_8
    aflux1(:,:,:,:,g) = 0.0_8
    sflux(:,:,g) = 0.0_8
    psii(:,:,g)= 0.0_8
    psij(:,:,g)= 0.0_8
  end do
  mu = 0.33_8
  eta = 0.66_8
  w = 0.4_8
  v = 0.1_8
  dx = 1.0_8 / nx
  dy = 1.0_8 / ny

  ! Allocate timers
  allocate(time(ntimes))

  if (rank.EQ.0) then
    write(*,'(a)') "MEGA-SWEEP!"
    write(*,*)
    write(*,'(a)') "Input"
    write(*,'(1x,a,i0,1x,a,1x,i0)') "Mesh size:          ", nx, "x", ny
    write(*,'(1x,a,i0)')      "Angles:             ", nang
    write(*,'(1x,a,i0)')      "Groups:             ", ng
    write(*,'(1x,a,i0)')      "Chunk:              ", chunk
    write(*,'(1x,a,i0)')      "Num. times:         ", ntimes
    write(*,*)
    write(*,'(a)') "Runtime info"
    write(*,'(1x,a,i0)')      "Num. procs:         ", nprocs
    write(*,'(1x,a,i0)')      "Num. threads/proc:  ", nthreads
    if (ydecomp) then
      write(*,'(1x,a)')       "Decomposing y-axis"
      if (mod(ny,nprocs) .eq. 0) then
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a)') "Sub-domain size:    ", nx, "x", lny, "(all ranks)"
      else
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", nx, "x", lny, "(", mod(ny,nprocs), "ranks)"
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", nx, "x", ny/nprocs, "(", nprocs-mod(ny,nprocs), "ranks)"
      end if
    else
      write(*,'(1x,a)')       "Decomposing x-axis"
      if (mod(nx,nprocs) .eq. 0) then
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a)') "Sub-domain size:    ", lnx, "x", lny, "(all ranks)"
      else
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", lnx, "x", ny, "(", mod(nx,nprocs), "ranks)"
        write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", nx/nprocs, "x", ny, "(", nprocs-mod(nx,nprocs), "ranks)"
      end if
    end if
    write(*,'(1x,a,f12.1)')   "Flux size (MB):     ", 8.0_8*(nang*nx*ny*nsweeps*ng)/2.0_8**20
    write(*,'(1x,a,f12.1)')   "Flux size/rank (MB):", 8.0_8*(nang*lnx*lny*nsweeps*ng)/2.0_8**20
    write(*,*)
  end if


  recv_time = 0.0_8
  wait_time = 0.0_8

  start_time = MPI_Wtime()
  do t = 1, ntimes

    timer = MPI_Wtime()

    ! Don't zero scalar flux, as it goes negative - OK for this benchmark with fake numbers

    if (ydecomp) then
      call sweeper_y(rank,lrank,rrank,           &
                   nang,nx,lny,ng,nsweeps,chunk, &
                   aflux0,aflux1,sflux,          &
                   psii,psij,                    &
                   mu,eta,w,v,dx,dy,buf)
    else
      call sweeper(rank,lrank,rrank,             &
                   nang,lnx,ny,ng,nsweeps,chunk, &
                   aflux0,aflux1,sflux,          &
                   psii,psij,                    &
                   mu,eta,w,v,dx,dy,buf)
    end if

    ! Swap pointers
    aflux_ptr => aflux0
    aflux0 => aflux1
    aflux1 => aflux_ptr

    time(t) = MPI_Wtime() - timer

  end do
  end_time = MPI_Wtime()

  total_time = end_time - start_time

  ! Collate slution
  call population(lnx,lny,ng,sflux,dx,dy,pop,total_pop)
  if (rank .eq. 0) then
    write(*,"(a)") "Population"
    write(*,"(1x,a,e23.16)") "Total:                   ", total_pop
    write(*,*)
  end if

  ! Model data movement
  moved = 8.0 * 1.0E-6 * (            &
          1.0*nang*nx*ny*ng*nsweeps + & ! read aflux0
          1.0*nang*nx*ny*ng*nsweeps + & ! write aflux1
          nang + nang +               & ! read mu and eta
          2.0*nang*ny*ng +            & ! read and write psii
          2.0*nang*nx*ng +            & ! read and write psij
          2.0*nx*ny*ng)                 ! read and write sflux

  lmoved = 8.0 * 1.0E-6 * (            &
          1.0*nang*lnx*lny*ng*nsweeps + & ! read aflux0
          1.0*nang*lnx*lny*ng*nsweeps + & ! write aflux1
          nang + nang +               & ! read mu and eta
          2.0*nang*lny*ng +            & ! read and write psii
          2.0*nang*lnx*ng +            & ! read and write psij
          2.0*lnx*lny*ng)                 ! read and write sflux

  tmoved = 8.0 * 1.0E-6 * (            &
          1.0*nang*lnx*lny*ng/nthreads*nsweeps + & ! read aflux0
          1.0*nang*lnx*lny*ng/nthreads*nsweeps + & ! write aflux1
          nang + nang +               & ! read mu and eta
          2.0*nang*lny*ng/nthreads +            & ! read and write psii
          2.0*nang*lnx*ng/nthreads +            & ! read and write psij
          2.0*lnx*lny*ng/nthreads)                 ! read and write sflux

  if (rank.EQ.0) then
    write(*,"(a)")   "Summary"
    write(*,"(1x,a,f12.9)") "Fastest iteration (s):   ", minval(time(2:))
    write(*,"(1x,a,f12.9)") "Slowest iteration (s)    ", maxval(time(2:))
    write(*,*)
    write(*,"(1x,a,f15.9,f5.1,a)") "Time in MPI_Recv (s):    ", recv_time, 100.0_8*recv_time/total_time, "%"
    write(*,"(1x,a,f15.9,f5.1,a)") "Time in MPI_Wait (s):    ", wait_time, 100.0_8*wait_time/total_time, "%"
    write(*,"(1x,a,f15.9,f5.1,a)") "Compute time (s):        ", &
      sum(time)-recv_time-wait_time, 100.0_8*(sum(time)-recv_time-wait_time)/total_time, "%"
    write(*,"(1x,a,f15.9,f5.1,a)") "Remaining time (s):      ", &
      total_time-sum(time), 100.0_8*(total_time-sum(time))/total_time, "%"
    write(*,*)
    write(*,"(1x,a,f15.9)") "Runtime (s):             ", total_time
    write(*,*)
    write(*,"(1x,a)")   "All ranks"
    write(*,"(2x,a,f12.2)") "Estimate moved (MB):     ", moved
    write(*,"(2x,a,f12.2)") "Best bandwidth (MB/s):   ", moved/minval(time(2:))
    write(*,"(2x,a,f12.2)") "Overall bandwidth (MB/s):", ntimes*moved/total_time
    write(*,*)
    write(*,"(1x,a)")   "Single rank"
    write(*,"(2x,a,f12.2)") "Best bandwidth (MB/s):   ", lmoved/minval(time(2:))
    write(*,"(2x,a,f12.2)") "Thread bandwidth (MB/s): ", tmoved/minval(time(2:))
    write(*,*)
    write(*,"(1x,a)")   "All ranks - cache bandwidth"
    write(*,"(2x,a,f12.2)") "Overall bandwidth (GB/s):", 1.0E-9*8*10*nang*nx*ny*ng*nsweeps*ntimes/total_time
  end if

  ! Free data
  deallocate(aflux0, aflux1)
  deallocate(sflux)
  deallocate(mu, eta)
  deallocate(psii, psij)
  deallocate(w)
  deallocate(pop)
  deallocate(time)

  call comms_finalize

end program

! Reduce the scalar flux
subroutine population(nx, ny, ng, sflux, dx, dy, pop, total)

  use comms

  implicit none

  integer :: nx, ny, ng
  real(kind=8) :: sflux(nx,ny,ng)
  real(kind=8) :: dx, dy
  real(kind=8) :: pop(ng)
  real(kind=8) :: total

  real(kind=8) :: tmp
  integer :: g

  do g = 1, ng
    tmp = sum(sflux(:,:,g)) * dx * dy
    call reduce(tmp, pop(g))
    pop(g) = pop(g) / real(ng - g + 1, 8)
  end do

  total = sum(pop)

end subroutine population

subroutine parse_args(rank,nang,nx,ny,ng,chunk,ntimes,ydecomp)

  implicit none

  integer, intent(in)   :: rank
  integer, intent(inout) :: nang, nx, ny, ng, chunk, ntimes
  logical, intent(inout) :: ydecomp

  character(len=32) :: arg

  integer :: i = 1

  do while (i <= command_argument_count())
    call getarg(i, arg)
    IF (arg .eq. "--nang") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) nang
    else if (arg .eq. "--ng") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) ng
    ELSE IF (arg .eq. "--nx") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) nx
    else if (arg .eq. "--ny") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) ny
    else if (arg .eq. "--chunk") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) chunk
    else if (arg .eq. "--ntimes") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) ntimes
    else if (arg .eq. "--ydecomp") then
      ydecomp = .true.
    else if (arg .eq. "--help") then
      if (rank .eq. 0) then
        write(*, *) "--nang   n  Set number of angles"
        write(*, *) "--ng     n  Set number of groups"
        write(*, *) "--nx     n  Set number of cells in x dimension"
        write(*, *) "--ny     n  Set number of cells in y dimension"
        write(*, *) "--chunk  n  Set y-dimension chunk size"
        write(*, *) "--ntimes n  Run the benchmark n times"
        write(*, *) "--ydecomp   Decompose in y instead of x (makes chunk over x-direction)"
        write(*, *)
        write(*, *) "Default sizes"
        write(*, '(2x,a,i0)') "nang: ", nang
        write(*, '(2x,a,i0)') "ng:   ", ng
        write(*, '(2x,a,i0)') "nx:   ", nx
        write(*, '(2x,a,i0)') "ny:   ", ny
        write(*, '(2x,a,i0)') "chunk:", chunk
      end if
      stop
    else
      if (rank .eq. 0) then
        write(*, *) "Unrecognised argument ", arg
        write(*, *)
      end if
      stop
    end if
    i = i + 1
  end do

end subroutine

