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


! Sweep kernel
subroutine sweeper(rank,lrank,rrank,            &
                   nang,nx,ny,ng,nsweeps,chunk, &
                   aflux0,aflux1,sflux,         &
                   psii,psij,                   &
                   mu,eta,                      &
                   w,v,dx,dy,                   &
                   buf,                         &
                   sweep_time,recv_time,send_time)

  use comms
  use MPI, only: MPI_Wtime

  implicit none

  integer :: rank, lrank, rrank
  integer :: nang, nx, ny, ng, nsweeps, chunk
  real(kind=8) :: aflux0(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: aflux1(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: sflux(nx,ny,ng)
  real(kind=8) :: psii(nang,chunk,ng)
  real(kind=8) :: psij(nang,nx,ng)
  real(kind=8) :: mu(nang)
  real(kind=8) :: eta(nang)
  real(kind=8) :: w(nang)
  real(kind=8) :: v
  real(kind=8) :: dx, dy
  real(kind=8) :: buf(nang,chunk,ng)
  real(kind=8) :: sweep_time(nsweeps)
  real(kind=8) :: recv_time, send_time

  integer :: a, i, j, g, c, cj, sweep
  integer :: istep, jstep ! Spatial step direction
  integer :: xmin, xmax   ! x-dimension loop bounds
  integer :: ymin, ymax   ! y-dimension (chunk) loop bounds
  integer :: cmin, cmax   ! Chunk loop bounds
  integer :: nchunks
  real(kind=8) :: psi
  real(kind=8) :: sweep_start, sweep_end
  real(kind=8) :: recv_start, recv_end
  real(kind=8) :: send_start, send_end

  ! Calculate number of chunks in y-dimension
  nchunks = ny / chunk

  do sweep = 1, nsweeps
    sweep_start = MPI_Wtime()

    ! Set sweep directions
    select case (sweep)
      case (1)
        istep = -1
        xmin = nx
        xmax = 1
        jstep = -1
        ymin = chunk
        ymax = 1
        cmin = nchunks
        cmax = 1
      case (2)
        istep = -1
        xmin = nx
        xmax = 1
        jstep = 1
        ymin = 1
        ymax = chunk
        cmin = 1
        cmax = nchunks
      case (3)
        istep = 1
        xmin = 1
        xmax = nx
        jstep = -1
        ymin = chunk
        ymax = 1
        cmin = nchunks
        cmax = 1
      case (4)
        istep = 1
        xmin = 1
        xmax = nx
        jstep = 1
        ymin = 1
        ymax = chunk
        cmin = 1
        cmax = nchunks
    end select

    ! Zero boundary data every sweep
    psii = 0.0_8
    psij = 0.0_8

    do c = cmin, cmax, jstep ! Loop over chunks

      ! Recv y boundary data for chunk
      recv_start = MPI_Wtime()
      psii = 0.0_8
      if (istep .eq. 1) then
        call recv(psii, nang*chunk*ng, lrank)
      else
        call recv(psii, nang*chunk*ng, rrank)
      end if
      recv_end = MPI_Wtime()
      recv_time = recv_time + recv_end - recv_start

      !$omp parallel do private(cj,j,i,a,psi)
      do g = 1, ng                 ! Loop over energy groups
        do cj = ymin, ymax, jstep  ! Loop over cells in chunk (y-dimension)
          ! Calculate y index with respect to ny
          j = (c-1)*chunk + cj
          do i = xmin, xmax, istep ! Loop over x-dimension
!dir$ vector nontemporal(aflux1)
            do a = 1, nang         ! Loop over angles
              ! Calculate angular flux
              psi = (mu(a)*psii(a,cj,g) + eta(a)*psij(a,i,g) + v*aflux0(a,i,j,sweep,g)) &
                    / (0.07_8 + 2.0_8*mu(a)/dx + 2.0_8*eta(a)/dy + v)

              ! Outgoing diamond difference
              psii(a,cj,g) = 2.0_8*psi - psii(a,cj,g)
              psij(a,i,g) = 2.0_8*psi - psij(a,i,g)
              aflux1(a,i,j,sweep,g) = 2.0_8*psi - aflux0(a,i,j,sweep,g)
  
              ! Reduction
              sflux(i,j,g) = sflux(i,j,g) + psi*w(a)

            end do ! angle loop
          end do ! x loop
        end do ! y chunk loop
      end do ! group loop
      !$omp end parallel do

      ! Send y boundary data for chunk
      ! NB non-blocking so need to buffer psii, making sure previous send has finished
      send_start = MPI_Wtime()
      call wait_on_sends
      buf = psii
      if (istep .eq. 1) then
        call send(buf, nang*chunk*ng, rrank)
      else
        call send(buf, nang*chunk*ng, lrank)
      end if
      send_end = MPI_Wtime()
      send_time = send_time + send_end - send_start

    end do ! chunk loop

    ! Total time in each sweep direction
    sweep_end = MPI_Wtime()
    sweep_time(sweep) = sweep_time(sweep) + sweep_end - sweep_start

  end do ! sweep loop

end subroutine sweeper

