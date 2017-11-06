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
  integer :: lr_send_request, ud_send_request

  ! Timer
  real(kind=8) :: recv_time, wait_time

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
    lr_send_request = MPI_REQUEST_NULL
    ud_send_request = MPI_REQUEST_NULL
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
  subroutine recv(lr_array,lr_num,lr_from,ud_array,ud_num,ud_from)

    real(kind=8) :: lr_array(:,:,:,:)
    real(kind=8) :: ud_array(:,:,:,:)
    integer, intent(in) :: lr_num, lr_from
    integer, intent(in) :: ud_num, ud_from
    integer :: err
    real(kind=8) :: time

    time = MPI_Wtime()
    call MPI_Recv(lr_array, lr_num, MPI_REAL8, lr_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    call MPI_Recv(ud_array, ud_num, MPI_REAL8, ud_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    recv_time = recv_time + MPI_Wtime() - time

  end subroutine recv

  ! Wait on previous send
  subroutine wait_on_sends

    integer :: err
    real(kind=8) :: time

    time = MPI_Wtime()
    call MPI_Wait(lr_send_request, MPI_STATUS_IGNORE, err)
    call MPI_Wait(ud_send_request, MPI_STATUS_IGNORE, err)
    wait_time = wait_time + MPI_Wtime() - time

  end subroutine wait_on_sends

  ! Send 4D arrays
  subroutine send(lr_array,lr_num,lr_to,ud_array,ud_num,ud_to)

    real(kind=8) :: lr_array(:,:,:,:)
    real(kind=8) :: ud_array(:,:,:,:)
    integer, intent(in) :: lr_num, lr_to
    integer, intent(in) :: ud_num, ud_to
    integer :: err

    ! wait_on_sends must have been previsouly called
    ! Array must be safe to send asynchronously, and not reused before the current Isend operation

    call MPI_Isend(lr_array, lr_num, MPI_REAL8, lr_to, 0, MPI_COMM_WORLD, lr_send_request, err)
    call MPI_Isend(ud_array, ud_num, MPI_REAL8, ud_to, 0, MPI_COMM_WORLD, ud_send_request, err)

  end subroutine send

  ! Reduce an array to single value
  subroutine reduce(val, res)

    real(kind=8) :: val, res
    integer :: err

    call MPI_Reduce(val, res, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, err)

  end subroutine

end module comms3d

