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

module comms3d

  use mpi

  implicit none

  ! MPI Send request, one for up/down and one for left/right
  integer :: y_send_request, z_send_request

  ! Timer
  real(kind=8) :: mpi_recv_time, mpi_wait_time

contains

  ! Init MPI
  subroutine comms_init

    integer :: mpi_thread_level
    integer :: err

    call MPI_Init_thread(MPI_THREAD_FUNNELED, mpi_thread_level, err)
    if (mpi_thread_level.LT.MPI_THREAD_FUNNELED) then
      print *, "Cannot provide MPI thread level"
      stop
    end if

    ! Set request to NULL
    y_send_request = MPI_REQUEST_NULL
    z_send_request = MPI_REQUEST_NULL
  end subroutine comms_init

  ! Finalize MPI
  subroutine comms_finalize

    integer :: err

    call MPI_Finalize(err)

  end subroutine comms_finalize

  ! Get rank number
  subroutine comms_rank(rank)

    integer, intent(out) :: rank
    integer :: err

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

  end subroutine comms_rank

  ! Get communicator size
  subroutine comms_size(nprocs)

    integer, intent(out) :: nprocs
    integer :: err

    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

  end subroutine comms_size

  ! foobar
  subroutine decompose(npey, npez, nprocs)

    integer, intent(out) :: npey, npez
    integer, intent(in) :: nprocs
    integer :: dims(2)
    integer :: err

    dims = 0
    call MPI_Dims_create(nprocs, 2, dims, err)
    npey = dims(1)
    npez = dims(2)

  end subroutine

  ! Receive 4D arrays
  subroutine recv(array,num,from)

    real(kind=8) :: array(:,:,:,:)
    integer, intent(in) :: num, from
    integer :: err
    real(kind=8) :: time

    time = MPI_Wtime()
    call MPI_Recv(array, num, MPI_REAL8, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    mpi_recv_time = mpi_recv_time + MPI_Wtime() - time

  end subroutine recv

  ! Wait on previous send
  subroutine wait_on_sends

    integer :: err
    real(kind=8) :: time

    time = MPI_Wtime()
    call MPI_Wait(y_send_request, MPI_STATUS_IGNORE, err)
    call MPI_Wait(z_send_request, MPI_STATUS_IGNORE, err)
    mpi_wait_time = mpi_wait_time + MPI_Wtime() - time

  end subroutine wait_on_sends

  ! Send 4D arrays
  subroutine ysend(array,num,to)

    real(kind=8) :: array(:,:,:,:)
    integer, intent(in) :: num, to
    integer :: err

    ! wait_on_sends must have been previsouly called
    ! Array must be safe to send asynchronously, and not reused before the current Isend operation

    call MPI_Isend(array, num, MPI_REAL8, to, 0, MPI_COMM_WORLD, y_send_request, err)

  end subroutine ysend

  ! Send 4D arrays
  subroutine zsend(array,num,to)

    real(kind=8) :: array(:,:,:,:)
    integer, intent(in) :: num, to
    integer :: err

    ! wait_on_sends must have been previsouly called
    ! Array must be safe to send asynchronously, and not reused before the current Isend operation

    call MPI_Isend(array, num, MPI_REAL8, to, 0, MPI_COMM_WORLD, z_send_request, err)

  end subroutine zsend

  ! Reduce an array to single value
  subroutine reduce(val, res)

    real(kind=8) :: val, res
    integer :: err

    call MPI_Reduce(val, res, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, err)

  end subroutine

end module comms3d

