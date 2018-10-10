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
program megasweep3d

  ! MEGA-SWEEP adds a KBA sweep along with the MPI communication to the MEGA-STREAM kernel
  ! 3D version

  use comms3d
  use omp_lib

  implicit none

  ! MPI variables
  integer :: rank, nprocs
  integer :: yrank, zrank
  integer :: npey, npez ! ranks in y and z dimensions
  integer :: ydown_rank, yup_rank, zdown_rank, zup_rank ! neighbour ranks
  real(kind=8), dimension(:,:,:,:), allocatable :: y_buf, z_buf ! comms buffers

  ! OpenMP variables
  integer :: nthreads

  ! Problem sizes
  integer :: nang, ng      ! angles and groups
  integer :: nx, ny, nz    ! global mesh size
  integer :: lnx, lny, lnz ! local mesh size
  integer :: chunk         ! chunk size
  integer :: nsweeps       ! sweep direction
  integer :: ntimes        ! number of times

  logical :: profile ! When true, turns of solution checking

  ! Arrays
  real(kind=8), dimension(:,:,:,:,:,:), pointer :: aflux0, aflux1   ! angular flux
  real(kind=8), dimension(:,:,:,:,:,:), pointer :: aflux_ptr        ! for pointer swap
  real(kind=8), dimension(:,:,:,:), allocatable :: sflux            ! scalar flux
  real(kind=8), dimension(:), allocatable :: mu, eta, xi            ! angular cosines
  real(kind=8), dimension(:,:,:,:), allocatable :: psii, psij, psik ! edge angular flux
  real(kind=8), dimension(:), allocatable :: w                      ! scalar flux weights
  real(kind=8) :: v                                                 ! time constant
  real(kind=8) :: dx, dy, dz                                        ! cell size
  real(kind=8), dimension(:), allocatable :: pop                    ! scalar flux reduction
  real(kind=8) :: total_pop

  ! Timers
  real(kind=8) :: start_time, end_time
  real(kind=8) :: total_time
  real(kind=8) :: timer
  real(kind=8), dimension(:), allocatable :: time
  real(kind=8), dimension(:), allocatable :: sweep_time
  real(kind=8) :: recv_time, send_time

  ! Local variables
  integer :: t, g, s
  real(kind=8) :: moved ! model of data movement

  call comms_init

  call comms_rank(rank)
  call comms_size(nprocs)

  !$omp parallel
    nthreads = omp_get_num_threads()
  !$omp end parallel

  ! Set default problem sizes
  nx = 128
  ny = 128
  nz = 128
  ng = 16
  nang = 16
  nsweeps = 8
  chunk = 1
  ntimes = 500

  profile = .false.

  ! Read in command line arguments
  call parse_args(rank,nang,nx,ny,nz,ng,chunk,ntimes,profile)

  ! Check nx can be split into even number of chunks
  if (mod(nx,chunk) .ne. 0) then
    if (rank .eq. 0) then
      print *, "Chunk size must divide nx"
    end if
    stop
  end if

  ! Decompose in y and z-dimension
  call decompose(npey, npez, nprocs)

  if (npey .gt. ny) then
    if (rank .eq. 0) then
      print *, "Too many processors to decompose ny"
    end if
    stop
  end if

  if (npez .gt. nz) then
    if (rank .eq. 0) then
      print *, "Too many processors to decompose nz"
    end if
    stop
  end if

  lnx = nx
  lny = ny / npey
  lnz = nz / npez

  ! Set ranks
  yrank = MOD(rank,npey)
  zrank = rank/npey

  ! Set neighbour ranks
  if (yrank .ne. 0) then
    ydown_rank = (yrank-1) + zrank*npey
  else
    ydown_rank = MPI_PROC_NULL
  end if
  if (yrank .ne. npey-1) then
    yup_rank = (yrank+1) + zrank*npey
  else
    yup_rank = MPI_PROC_NULL
  end if
  if (zrank .ne. 0) then
    zdown_rank = yrank + (zrank-1)*npey
  else
    zdown_rank = MPI_PROC_NULL
  end if
  if (zrank .ne. npez-1) then
    zup_rank = yrank + (zrank+1)*npey
  else
    zup_rank = MPI_PROC_NULL
  end if

  ! Share remainder cells for uneven decomposition
  ! Allows for flexible process counts
  if (mod(ny,npey) .ne. 0) then
    if (yrank .lt. mod(ny,npey)) then
      lny = lny + 1
    end if
  end if
  if (mod(nz,npez) .ne. 0) then
    if (zrank .lt. mod(nz,npez)) then
      lnz = lnz + 1
    end if
  end if

  if (rank.EQ.0) then
    write(*,'(a)') "MEGA-SWEEP!"
    write(*,*)
    write(*,'(a)') "Input"
    write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,1x,i0)') "Mesh size:          ", nx, "x", ny, "x", nz
    write(*,'(1x,a,i0)')      "Angles:             ", nang
    write(*,'(1x,a,i0)')      "Groups:             ", ng
    write(*,'(1x,a,i0)')      "Chunk:              ", chunk
    write(*,'(1x,a,i0)')      "Num. times:         ", ntimes
    write(*,*)
    write(*,'(a)') "Runtime info"
    write(*,'(1x,a,i0)')      "Num. procs:         ", nprocs
    write(*,'(1x,a,i0)')      "Num. threads/proc:  ", nthreads
    write(*,'(1x,a)')       "Decomposing x-axis"
    write(*,'(1x,a,i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') "Sub-domain size:    ", lnx, "x", lny, "x", lnz, "(master)"
    if (mod(ny,npey) .ne. 0 .or. mod(nz,npez) .ne. 0) then
      write(*,'(1x,a)') "Note: Uneven decomposition"
    end if
    write(*,'(1x,a,f12.1)')   "Flux size (MB):     ", (8.0_8*nang*nx*ny*nz*nsweeps*ng)/2.0_8**20
    write(*,'(1x,a,f12.1)')   "Flux size/rank (MB):", (8.0_8*nang*lnx*lny*lnz*nsweeps*ng)/2.0_8**20
    write(*,*)
  end if

  ! Allocate data
  allocate(aflux0(nang,lnx,lny,lnz,nsweeps,ng))
  allocate(aflux1(nang,lnx,lny,lnz,nsweeps,ng))
  allocate(sflux(lnx,lny,lnz,ng))
  allocate(mu(nang))
  allocate(eta(nang))
  allocate(xi(nang))
  allocate(psii(nang,lny,lnz,ng))
  allocate(psij(nang,chunk,lnz,ng))
  allocate(psik(nang,chunk,lny,ng))
  allocate(w(nang))
  allocate(pop(ng))
  allocate(y_buf(nang,chunk,lnz,ng))
  allocate(z_buf(nang,chunk,lny,ng))

  ! Allocate timers
  allocate(time(ntimes))
  allocate(sweep_time(nsweeps))

  ! Initilise data
  !$omp parallel do
  do g = 1, ng
    aflux0(:,:,:,:,:,g) = 1.0_8
    aflux1(:,:,:,:,:,g) = 0.0_8
    sflux(:,:,:,g) = 0.0_8
    psii(:,:,:,g)= 0.0_8
    psij(:,:,:,g)= 0.0_8
    psik(:,:,:,g)= 0.0_8
  end do
  mu = 0.33_8
  eta = 0.66_8
  xi = sqrt(1.0_8 - mu(1)*mu(1) - eta(1)*eta(1))
  w = 0.4_8
  v = 0.1_8
  dx = 1.0_8 / nx
  dy = 1.0_8 / ny
  dz = 1.0_8 / nz

  mpi_recv_time = 0.0_8
  mpi_wait_time = 0.0_8
  sweep_time = 0.0_8
  send_time = 0.0_8
  recv_time = 0.0_8

  start_time = MPI_Wtime()
  do t = 1, ntimes

    timer = MPI_Wtime()

    ! Don't zero scalar flux, as it goes negative - OK for this benchmark with fake numbers

    call sweeper3d(rank,yup_rank,ydown_rank,zup_rank,zdown_rank,      &
                   nang,lnx,lny,lnz,ng,nsweeps,chunk, &
                   aflux0,aflux1,sflux,               &
                   psii,psij,psik,                    &
                   mu,eta,xi,w,v,dx,dy,dz,            &
                   y_buf,z_buf,                       &
                   sweep_time,recv_time,send_time)


    ! Swap pointers
    aflux_ptr => aflux0
    aflux0 => aflux1
    aflux1 => aflux_ptr

    time(t) = MPI_Wtime() - timer

  end do
  end_time = MPI_Wtime()

  total_time = end_time - start_time

  if (.not.profile) then

    ! Check angular flux is non-zero
    ! If still zero, then sweep didn't touch a cell
    if (any(aflux0 .eq. 0.0_8)) then
      write(*,"(a,i0)") "Warning: angular flux 0 contains zero values on rank ", rank
      write(*,*)
    end if
    if (any(aflux1 .eq. 0.0_8)) then
      write(*,"(a,i0)") "Warning: angular flux 1 contains zero values on rank ", rank
      write(*,*)
    end if

    ! Collate solution
    call population(lnx,lny,lnz,ng,sflux,dx,dy,dz,pop,total_pop)
    if (rank .eq. 0) then
      write(*,"(a)") "Population"
      do g = 1, ng
        write(*,"(1x,i0,a,e23.16)") g, ":", pop(g)
      end do
      write(*,"(1x,a,e23.16)") "Total:                   ", total_pop
      write(*,*)
    end if

  end if


  ! Model data movement
  moved = 8.0 * 1.0E-6 * (            &
          1.0*nang*nx*ny*nz*ng*nsweeps + & ! read aflux0
          1.0*nang*nx*ny*nz*ng*nsweeps + & ! write aflux1
          nang + nang + nang +        & ! read mu and eta and xi
          2.0*nang*ny*nz*ng +            & ! read and write psii
          2.0*nang*nx*nz*ng +            & ! read and write psij
          2.0*nang*nx*ny*ng +            & ! read and write psik
          2.0*nx*ny*nz*ng)                 ! read and write sflux

  if (rank.EQ.0) then
    write(*,"(a)")   "Summary"
    write(*,"(1x,a,f12.9)") "Fastest iteration (s):   ", minval(time(2:))
    write(*,"(1x,a,f12.9)") "Slowest iteration (s)    ", maxval(time(2:))
    write(*,*)
    write(*,"(1x,a)") "Timings (s)"
    write(*,"(2x,a,f15.9)") "Runtime:         ", total_time
    write(*,"(2x,a,f15.9)") "Solve time:      ", sum(time)
    do s = 1, nsweeps
      write(*,"(3x,a,i0,a,f15.9,f5.1,a)") "Sweep ", s, ":        ", sweep_time(s), sweep_time(s)/total_time*100.0_8, "%"
    end do
    write(*,*)

    write(*,"(3x,a,f15.9,f5.1,a)") "Compute time:   ", &
      sum(time)-recv_time-send_time, 100.0_8*(sum(time)-recv_time-send_time)/total_time, "%"
    write(*,*)

    write(*,"(3x,a,f15.9,f5.1,a)") "Communication:  ", recv_time+send_time, (recv_time+send_time)/total_time*100.0_8, "%"
    write(*,"(4x,a,f15.9,f5.1,a)")  "Receives:      ", recv_time, recv_time/total_time*100.0_8, "%"
    write(*,"(5x,a,f15.9,f5.1,a)")   "MPI_Recv:     ", mpi_recv_time, 100.0_8*mpi_recv_time/total_time, "%"
    write(*,"(4x,a,f15.9,f5.1,a)")  "Sends:         ", send_time, send_time/total_time*100.0_8, "%"
    write(*,"(5x,a,f15.9,f5.1,a)")   "MPI_Wait:     ", mpi_wait_time, 100.0_8*mpi_wait_time/total_time, "%"
    write(*,*)

    write(*,"(2x,a,f15.9,f5.1,a)") "Remaining time:  ", &
      total_time-sum(time), 100.0_8*(total_time-sum(time))/total_time, "%"
    write(*,*)

    write(*,"(1x,a)")   "All ranks"
    write(*,"(2x,a,f12.2)") "Estimate moved (MB):     ", moved
    write(*,"(2x,a,f12.2)") "Best bandwidth (MB/s):   ", moved/minval(time(2:))
    write(*,"(2x,a,f12.2)") "Overall bandwidth (MB/s):", ntimes*moved/total_time
    write(*,*)
  end if

  ! Free data
  deallocate(aflux0, aflux1)
  deallocate(sflux)
  deallocate(mu, eta, xi)
  deallocate(psii, psij, psik)
  deallocate(w)
  deallocate(pop)
  deallocate(time)
  deallocate(sweep_time)

  call comms_finalize

end program

! Reduce the scalar flux
subroutine population(nx, ny, nz, ng, sflux, dx, dy, dz, pop, total)

  use comms3d

  implicit none

  integer :: nx, ny, nz, ng
  real(kind=8) :: sflux(nx,ny,nz,ng)
  real(kind=8) :: dx, dy, dz
  real(kind=8) :: pop(ng)
  real(kind=8) :: total

  real(kind=8) :: tmp
  integer :: g

  do g = 1, ng
    tmp = sum(sflux(:,:,:,g)) * dx * dy * dz
    call reduce(tmp, pop(g))
    pop(g) = pop(g) / real(ng - g + 1, 8)
  end do

  total = sum(pop)

end subroutine population

subroutine parse_args(rank,nang,nx,ny,nz,ng,chunk,ntimes,profile)

  implicit none

  integer, intent(in)   :: rank
  integer, intent(inout) :: nang, nx, ny, nz, ng, chunk, ntimes
  logical, intent(inout) :: profile

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
    else if (arg .eq. "--nz") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) nz
    else if (arg .eq. "--chunk") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) chunk
    else if (arg .eq. "--ntimes") then
      i = i + 1
      call getarg(i, arg)
      read(arg, *) ntimes
    else if (arg .eq. "--prof") then
      profile = .true.
    else if (arg .eq. "--help") then
      if (rank .eq. 0) then
        write(*, *) "--nang   n  Set number of angles"
        write(*, *) "--ng     n  Set number of groups"
        write(*, *) "--nx     n  Set number of cells in x dimension"
        write(*, *) "--ny     n  Set number of cells in y dimension"
        write(*, *) "--nz     n  Set number of cells in z dimension"
        write(*, *) "--chunk  n  Set x-dimension chunk size"
        write(*, *) "--ntimes n  Run the benchmark n times"
        write(*,*)  "--prof      Run in profiler mode: removes solution check"
        write(*, *)
        write(*, *) "Default sizes"
        write(*, '(2x,a,i0)') "nang: ", nang
        write(*, '(2x,a,i0)') "ng:   ", ng
        write(*, '(2x,a,i0)') "nx:   ", nx
        write(*, '(2x,a,i0)') "ny:   ", ny
        write(*, '(2x,a,i0)') "chunk:", chunk
        write(*,*)
        write(*, '(2x,a)') "MPI decomposition is in YZ dimensions."
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

