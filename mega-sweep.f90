! mega-stream within a KBA sweep
program megasweep

  use mpi

  implicit none

  ! MPI variables
  integer :: mpi_thread_level
  integer :: ierr
  integer :: rank, nprocs
  integer :: comm


  ! Problem sizes
  integer :: nang, ng ! angles and groups
  integer :: nx, ny   ! global mesh size
  integer :: lnx      ! local mesh size
  integer :: chunk    ! y chunk size
  integer :: nsweeps  ! sweep direction

  ! Arrays
  real(kind=8), dimension(:,:,:,:,:), pointer :: aflux0, aflux1 ! angular flux
  real(kind=8), dimension(:,:,:,:,:), pointer :: aflux_ptr      ! for pointer swap
  real(kind=8), dimension(:,:,:), allocatable :: sflux          ! scalar flux
  real(kind=8), dimension(:), allocatable :: mu, eta            ! angular cosines
  real(kind=8), dimension(:,:,:), allocatable :: psii, psij      ! edge angular flux

  call MPI_Init_thread(MPI_THREAD_FUNNELED, mpi_thread_level, ierr)
  if (mpi_thread_level.LT.MPI_THREAD_FUNNELED) then
    print *, "Cannot provide MPI thread level"
    stop
  end if

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nprocs, ierr)

  ! Set problem sizes
  nx = 128
  ny = 128
  ng = 16
  nang = 16
  nsweeps = 4
  chunk = 1

  ! Decompose in x-dimension
  if (mod(nx,nprocs).NE.0) then
    if (rank.EQ.0) then
      print *, "Number of ranks must divide nx"
    end if
    stop
  end if
  lnx = nx / nprocs

  ! Allocate data
  allocate(aflux0(nang,nx,ny,nsweeps,ng))
  allocate(aflux1(nang,nx,ny,nsweeps,ng))
  allocate(sflux(nx,ny,ng))
  allocate(mu(nang))
  allocate(eta(nang))
  allocate(psii(nang,chunk,ng))
  allocate(psij(nang,nx,ng))

  if (rank.EQ.0) then
    print *, "MEGA-SWEEP!"
    print *
    print *, "Num. procs:", nprocs
    print *, "Mesh size:", nx, ny
  end if

  call sweep(nang,lnx,ny,ng,nsweeps,chunk, &
              aflux0,aflux1,sflux,          &
              psii,psij,                    &
              mu,eta)

  ! Free data
  deallocate(aflux0, aflux1)
  deallocate(sflux)
  deallocate(mu, eta)
  deallocate(psii, psij)

  call MPI_Finalize(ierr)

end program

! Sweep kernel
subroutine sweep(nang,nx,ny,ng,nsweeps,chunk, &
                 aflux0,aflux1,sflux,         &
                 psii,psij,                   &
                 mu,eta)

  implicit none

  integer :: nang, nx, ny, ng, nsweeps, chunk
  real(kind=8) :: aflux0(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: aflux1(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: sflux(nx,ny,ng)
  real(kind=8) :: psii(nang,chunk,ng)
  real(kind=8) :: psij(nang,ny,ng)
  real(kind=8) :: mu(nang)
  real(kind=8) :: eta(nang)

  print *, shape(aflux0)

end subroutine sweep

