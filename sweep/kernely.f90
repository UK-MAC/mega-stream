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


! Sweep kernel, but mesh is decomposed in y-dimension
subroutine sweeper_y(rank,lrank,rrank,          &
                   nang,nx,ny,ng,nsweeps,chunk, &
                   aflux0,aflux1,sflux,         &
                   psii,psij,                   &
                   mu,eta,                      &
                   w,v)

  use comms

  implicit none

  integer :: rank, lrank, rrank
  integer :: nang, nx, ny, ng, nsweeps, chunk
  real(kind=8) :: aflux0(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: aflux1(nang,nx,ny,nsweeps,ng)
  real(kind=8) :: sflux(nx,ny,ng)
  real(kind=8) :: psii(nang,ny,ng)
  real(kind=8) :: psij(nang,chunk,ng)
  real(kind=8) :: mu(nang)
  real(kind=8) :: eta(nang)
  real(kind=8) :: w(nang)
  real(kind=8) :: v

  integer :: a, i, j, g, c, ci, sweep
  integer :: istep, jstep ! Spatial step direction
  integer :: xmin, xmax   ! x-dimension (chunk) loop bounds
  integer :: ymin, ymax   ! y-dimension loop bounds
  integer :: cmin, cmax   ! Chunk loop bounds
  integer :: nchunks
  real(kind=8) :: psi

  ! Calculate number of chunks in x-dimension
  nchunks = nx / chunk

  do sweep = 1, nsweeps
    ! Set sweep directions
    select case (sweep)
      case (1)
        istep = -1
        xmin = chunk
        xmax = 1
        cmin = nchunks
        cmax = 1
        jstep = -1
        ymin = ny
        ymax = 1
      case (2)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = -1
        ymin = ny
        ymax = 1
      case (3)
        istep = -1
        xmin = chunk
        xmax = 1
        cmin = nchunks
        cmax = 1
        jstep = 1
        ymin = 1
        ymax = ny
      case (4)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = 1
        ymin = 1
        ymax = ny
    end select

    ! Zero boundary data every sweep
    psii = 0.0_8
    psij = 0.0_8

    do c = cmin, cmax, istep ! Loop over chunks

      ! Recv x boundary data for chunk
      psij = 0.0_8
      if (jstep .eq. 1) then
        call recv(psij, nang*chunk*ng, rrank)
      else
        call recv(psij, nang*chunk*ng, lrank)
      end if

      do g = 1, ng                 ! Loop over energy groups
        do j = ymin, ymax, jstep  ! Loop over y-dimension
          do ci = xmin, xmax, istep ! Loop over cells in chunk (x-dimension)
            ! Calculate x index with respect to nx
            i = (c-1)*chunk + ci
!dir$ vector nontemporal(aflux1)
            do a = 1, nang         ! Loop over angles
              ! Calculate angular flux
              psi = mu(a)*psii(a,j,g) + eta(a)*psij(a,ci,g) + v*aflux0(a,i,j,sweep,g)

              ! Outgoing diamond difference
              psii(a,j,g) = 2.0_8*psi - psii(a,j,g)
              psij(a,ci,g) = 2.0_8*psi - psij(a,ci,g)
              aflux1(a,i,j,sweep,g) = 2.0_8*psi - aflux0(a,i,j,sweep,g)
  
              ! Reduction
              sflux(i,j,g) = sflux(i,j,g) + psi*w(a)

            end do ! angle loop
          end do ! x loop
        end do ! y chunk loop
      end do ! group loop

      ! Send y boundary data for chunk
      if (jstep .eq. 1) then
        call send(psij, nang*chunk*ng, lrank)
      else
        call send(psij, nang*chunk*ng, rrank)
      end if

    end do ! chunk loop
  end do ! sweep loop

end subroutine sweeper_y

