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

  implicit none

  ! MPI variables
  integer :: rank, nprocs
  integer :: lrank, rrank ! neighbour ranks


  ! Problem sizes
  integer :: nang, ng ! angles and groups
  integer :: nx, ny   ! global mesh size
  integer :: lnx      ! local mesh size
  integer :: chunk    ! y chunk size
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

  ! Timers
  real(kind=8) :: start_time, end_time
  real(kind=8) :: timer
  real(kind=8), dimension(:), allocatable :: time

  ! Local variables
  integer :: t
  real(kind=8) :: moved ! model of data movement

  call comms_init

  call comms_rank(rank)
  call comms_size(nprocs)

  ! Set default problem sizes
  nx = 128
  ny = 128
  ng = 16
  nang = 16
  nsweeps = 4
  chunk = 1
  ntimes = 500

  ! Read in command line arguments
  call parse_args(rank,nang,nx,ny,ng,chunk,ntimes)

  ! Check ny is split into even number of chunks
  if (mod(ny,chunk) .ne. 0) then
    if (rank .eq. 0) then
      print *, "Chunk size must divide ny"
    end if
    stop
  end if

  ! Decompose in x-dimension
  lnx = nx / nprocs

  ! Share remainder cells for uneven decomposition 
  ! Allows for flexible process counts
  if (mod(nx,nprocs) .ne. 0) then
    if (rank .lt. mod(nx,nprocs)) then
      lnx = lnx + 1
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
  allocate(aflux0(nang,lnx,ny,nsweeps,ng))
  allocate(aflux1(nang,lnx,ny,nsweeps,ng))
  allocate(sflux(lnx,ny,ng))
  allocate(mu(nang))
  allocate(eta(nang))
  allocate(psii(nang,chunk,ng))
  allocate(psij(nang,lnx,ng))
  allocate(w(nang))

  ! Initilise data
  aflux0 = 1.0_8
  aflux1 = 0.0_8
  sflux = 0.0_8
  mu = 0.33_8
  eta = 0.66_8
  psii = 0.0_8
  psij = 0.0_8
  w = 0.4_8
  v = 0.1_8

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
    if (mod(nx,nprocs) .eq. 0) then
      write(*,'(1x,a,i0,1x,a,1x,i0,1x,a)') "Sub-domain size:    ", lnx, "x", ny, "(all ranks)"
    else
      write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", lnx, "x", ny, "(", mod(nx,nprocs), "ranks)"
      write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,i0,1x,a)') "Sub-domain size:    ", nx/nprocs, "x", ny, "(", nprocs-mod(nx,nprocs), "ranks)"
    end if
    write(*,'(1x,a,f12.1)')   "Flux size (MB):     ", 8.0_8*(nang*nx*ny*nsweeps*ng)/2.0_8**20
    write(*,'(1x,a,f12.1)')   "Flux size/rank (MB):", 8.0_8*(nang*lnx*ny*nsweeps*ng)/2.0_8**20
    write(*,*)
  end if


  start_time = MPI_Wtime()
  do t = 1, ntimes

    timer = MPI_Wtime()

    call sweeper(rank,lrank,rrank,             &
                 nang,lnx,ny,ng,nsweeps,chunk, &
                 aflux0,aflux1,sflux,          &
                 psii,psij,                    &
                 mu,eta,w,v)

    ! Swap pointers
    aflux_ptr => aflux0
    aflux0 => aflux1
    aflux1 => aflux_ptr

    time(t) = MPI_Wtime() - timer

  end do
  end_time = MPI_Wtime()

  ! Model data movement
  moved = 8.0 * 1.0E-6 * (            &
          1.0*nang*nx*ny*ng*nsweeps + & ! read aflux0
          1.0*nang*nx*ny*ng*nsweeps + & ! write aflux1
          nang + nang +               & ! read mu and eta
          2.0*nang*ny*ng +            & ! read and write psii
          2.0*nang*nx*ng +            & ! read and write psij
          2.0*nx*ny*ng)                 ! read and write sflux


  if (rank.EQ.0) then
    write(*,"(a)")   "Summary"
    write(*,"(1x,a,f12.9)") "Fastest iteration (s):   ", minval(time(2:))
    write(*,"(1x,a,f12.9)") "Slowest iteration (s)    ", maxval(time(2:))
    write(*,"(1x,a,f12.9)") "Runtime (s):             ", end_time-start_time
    write(*,*)
    write(*,"(1x,a,f12.2)") "Estimate moved (MB):     ", moved
    write(*,*)
    write(*,"(1x,a,f12.2)") "Best bandwidth (MB/s):   ", moved/minval(time(2:))
    write(*,"(1x,a,f12.2)") "Overall bandwidth (MB/s):", ntimes*moved/(end_time-start_time)
    write(*,*)
  end if

  ! Free data
  deallocate(aflux0, aflux1)
  deallocate(sflux)
  deallocate(mu, eta)
  deallocate(psii, psij)
  deallocate(w)
  deallocate(time)

  call comms_finalize

end program

subroutine parse_args(rank,nang,nx,ny,ng,chunk,ntimes)

  implicit none

  integer, intent(in)   :: rank
  integer, intent(inout) :: nang, nx, ny, ng, chunk, ntimes

  character(len=32) :: arg

  integer :: i = 1

  do while (i <= iargc())
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
    else if (arg .eq. "--help") then
      if (rank .eq. 0) then
        write(*, *) "--nang   n  Set number of angles"
        write(*, *) "--ng     n  Set number of groups"
        write(*, *) "--nx     n  Set number of cells in x dimension"
        write(*, *) "--ny     n  Set number of cells in y dimension"
        write(*, *) "--chunk  n  Set y-dimension chunk size"
        write(*, *) "--ntimes n  Run the benchmark n times"
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

