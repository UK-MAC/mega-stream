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
                   w,v)

  use comms

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

  integer :: a, i, j, g, cj, sweep
  integer :: istep, jstep, xmin, xmax, ymin, ymax, cmin, cmax
  real(kind=8) :: psi

  do sweep = 1, nsweeps
    ! Set sweep directions
    select case (sweep)
      case (1)
        istep = -1
        xmin = nx
        xmax = 1
        jstep = -1
        ymin = ny
        ymax = 1
        cmin = chunk
        cmax = 1
      case (2)
        istep = 1
        xmin = 1
        xmax = nx
        jstep = -1
        ymin = ny
        ymax = 1
        cmin = chunk
        cmax = 1
      case (3)
        istep = -1
        xmin = nx
        xmax = 1
        jstep = 1
        ymin = 1
        ymax = ny
        cmin = 1
        cmax = chunk
      case (4)
        istep = 1
        xmin = 1
        xmax = nx
        jstep = 1
        ymin = 1
        ymax = ny
        cmin = 1
        cmax = chunk
    end select

    ! Zero boundary data every sweep
    psii = 0.0_8
    psij = 0.0_8

    do j = ymin, ymax, chunk*jstep

      ! Recv y boundary data for chunk
      if (istep .eq. 1) then
        call recv(psii, nang*chunk*ng, lrank)
      else
        call recv(psii, nang*chunk*ng, rrank)
      end if

      do g = 1, ng
        do cj = cmin, cmax, jstep
          do i = xmin, xmax, istep
!dir$ vector nontemporal(aflux1)
            do a = 1, nang
              ! Calculate angular flux
              psi = mu(a)*psii(a,cj,g) + eta(a)*psij(a,i,g) + v*aflux0(a,i,j+cj-1,sweep,g)

              ! Outgoing diamond difference
              psii(a,cj,g) = 2.0_8*psi - psii(a,cj,g)
              psij(a,i,g) = 2.0_8*psi - psij(a,i,g)
              aflux1(a,i,j+cj-1,sweep,g) = 2.0_8*psi - aflux0(a,i,j+cj-1,sweep,g)
  
              ! Reduction
              sflux(i,j+cj-1,g) = sflux(i,j+cj-1,g) + psi*w(a)

            end do ! angle loop
          end do ! x loop
        end do ! y chunk loop
      end do ! group loop

      ! Send y boundary data for chunk
      if (istep .eq. 1) then
        call send(psii, nang*chunk*ng, rrank)
      else
        call send(psii, nang*chunk*ng, lrank)
      end if

    end do ! chunk loop
  end do ! sweep loop

end subroutine sweeper

