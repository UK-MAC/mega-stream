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
subroutine sweeper3d(rank,yup_rank,ydown_rank,zup_rank,zdown_rank,      &
                     nang,nx,ny,nz,ng,nsweeps,chunk, &
                     aflux0,aflux1,sflux,               &
                     psii,psij,psik,                    &
                     mu,eta,xi,w,v,dx,dy,dz,            &
                     y_buf,z_buf,                       &
                     sweep_time,recv_time,send_time)

  use comms3d
  use MPI, only: MPI_Wtime

  implicit none

  integer :: rank, yup_rank, ydown_rank, zup_rank, zdown_rank
  integer :: nang, nx, ny, nz, ng, nsweeps, chunk
  real(kind=8) :: aflux0(nang,nx,ny,nz,nsweeps,ng)
  real(kind=8) :: aflux1(nang,nx,ny,nz,nsweeps,ng)
  real(kind=8) :: sflux(nx,ny,nz,ng)
  real(kind=8) :: psii(nang,ny,nz,ng)
  real(kind=8) :: psij(nang,chunk,nz,ng)
  real(kind=8) :: psik(nang,chunk,ny,ng)
  real(kind=8) :: mu(nang)
  real(kind=8) :: eta(nang)
  real(kind=8) :: xi(nang)
  real(kind=8) :: w(nang)
  real(kind=8) :: v
  real(kind=8) :: dx, dy, dz
  real(kind=8) :: y_buf(nang,chunk,nz,ng)
  real(kind=8) :: z_buf(nang,chunk,ny,ng)
  real(kind=8) :: sweep_time(nsweeps)
  real(kind=8) :: recv_time, send_time

  integer :: a, i, j, k, g, c, ci, sweep
  integer :: istep, jstep, kstep ! Spatial step direction
  integer :: xmin, xmax   ! x-dimension loop bounds
  integer :: ymin, ymax   ! y-dimension loop bounds
  integer :: zmin, zmax   ! z-dimension loop bounds
  integer :: cmin, cmax   ! Chunk loop bounds
  integer :: nchunks
  real(kind=8) :: psi
  real(kind=8) :: sweep_start, sweep_end
  real(kind=8) :: invdx, invdy, invdz
  real(kind=8) :: recv_start, recv_end
  real(kind=8) :: send_start, send_end

  ! Calculate number of chunks in x-dimension
  nchunks = nx / chunk

  ! Invert spatial width to remove divisions
  invdx = 2.0_8 / dx
  invdy = 2.0_8 / dy
  invdz = 2.0_8 / dz

  do sweep = 1, nsweeps
    sweep_start = MPI_Wtime()

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
        kstep = -1
        zmin = nz
        zmax = 1
      case (2)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = -1
        ymin = ny
        ymax = 1
        kstep = -1
        zmin = nz
        zmax = 1
      case (3)
        istep = -1
        xmin = chunk
        xmax = 1
        cmin = nchunks
        cmax = 1
        jstep = 1
        ymin = 1
        ymax = ny
        kstep = -1
        zmin = nz
        zmax = 1
      case (4)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = 1
        ymin = 1
        ymax = ny
        kstep = -1
        zmin = nz
        zmax = 1
      case (5)
        istep = -1
        xmin = chunk
        xmax = 1
        cmin = nchunks
        cmax = 1
        jstep = -1
        ymin = ny
        ymax = 1
        kstep = 1
        zmin = 1
        zmax = nz
      case (6)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = -1
        ymin = ny
        ymax = 1
        kstep = 1
        zmin = 1
        zmax = nz
      case (7)
        istep = -1
        xmin = chunk
        xmax = 1
        cmin = nchunks
        cmax = 1
        jstep = 1
        ymin = 1
        ymax = ny
        kstep = 1
        zmin = 1
        zmax = nz
      case (8)
        istep = 1
        xmin = 1
        xmax = chunk
        cmin = 1
        cmax = nchunks
        jstep = 1
        ymin = 1
        ymax = ny
        kstep = 1
        zmin = 1
        zmax = nz
    end select

    ! Zero boundary data every sweep
    !$omp parallel do
    do g = 1, ng
      psii(:,:,:,g) = 0.0_8
      psij(:,:,:,g) = 0.0_8
      psik(:,:,:,g) = 0.0_8
    end do
    !$omp end parallel do

    do c = cmin, cmax, istep ! Loop over chunks

      ! Recv boundary data for chunk
      recv_start = MPI_Wtime()

      !$omp parallel do
      do g = 1, ng
        psij(:,:,:,g) = 0.0_8
        psik(:,:,:,g) = 0.0_8
      end do
      !$omp end parallel do

      if (jstep .eq. 1) then
        call recv(psij, nang*chunk*nz*ng, ydown_rank)
      else
        call recv(psij, nang*chunk*nz*ng, yup_rank)
      end if
      if (kstep .eq. 1) then
        call recv(psik, nang*chunk*ny*ng, zdown_rank)
      else
        call recv(psik, nang*chunk*ny*ng, zup_rank)
      end if
      recv_end = MPI_Wtime()
      recv_time = recv_time + recv_end - recv_start

      !$omp parallel do private(k,j,ci,i,a,psi)
      do g = 1, ng                 ! Loop over energy groups
        do k = zmin, zmax, kstep   ! Loop over z-dimension
          do j = ymin, ymax, jstep ! Loop over y-dimension
            do ci = xmin, xmax, istep  ! Loop over cells in chunk (x-dimension)
              ! Calculate x index with respect to nx
              i = (c-1)*chunk + ci
              !dir$ vector nontemporal(aflux1)
              do a = 1, nang         ! Loop over angles
                ! Calculate angular flux
                psi = (mu(a)*psii(a,j,k,g) + eta(a)*psij(a,ci,k,g) + xi(a)*psik(a,ci,j,g) + v*aflux0(a,i,j,k,sweep,g)) &
                      / (0.07_8 + mu(a)*invdx + eta(a)*invdy + xi(a)*invdz + v)

                ! Outgoing diamond difference
                psii(a,j,k,g) = 2.0_8*psi - psii(a,j,k,g)
                psij(a,ci,k,g) = 2.0_8*psi - psij(a,ci,k,g)
                psik(a,ci,j,g) = 2.0_8*psi - psik(a,ci,j,g)
                aflux1(a,i,j,k,sweep,g) = 2.0_8*psi - aflux0(a,i,j,k,sweep,g)
  
                ! Reduction
                sflux(i,j,k,g) = sflux(i,j,k,g) + psi*w(a)

              end do ! angle loop
            end do ! x chunk loop
          end do ! y loop
        end do ! z loop
      end do ! group loop
      !$omp end parallel do

      ! Send y boundary data for chunk
      ! NB non-blocking so need to buffer psii, making sure previous send has finished
      send_start = MPI_Wtime()
      call wait_on_sends

      !$omp parallel do
      do g = 1, ng
        y_buf(:,:,:,g) = psij(:,:,:,g)
        z_buf(:,:,:,g) = psik(:,:,:,g)
      end do
      !$omp end parallel do

      if (jstep .eq. 1) then
        call ysend(y_buf, nang*chunk*nz*ng, yup_rank)
      else
        call ysend(y_buf, nang*chunk*nz*ng, ydown_rank)
      end if
      if (kstep .eq. 1) then
        call zsend(z_buf, nang*chunk*ny*ng, zup_rank)
      else
        call zsend(z_buf, nang*chunk*ny*ng, zdown_rank)
      end if
      send_end = MPI_Wtime()
      send_time = send_time + send_end - send_start

    end do ! chunk loop

    ! Total time in each sweep direction
    sweep_end = MPI_Wtime()
    sweep_time(sweep) = sweep_time(sweep) + sweep_end - sweep_start

  end do ! sweep loop

end subroutine sweeper3d

