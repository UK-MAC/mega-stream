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

! This is the Fortran interface to envoke PAPI via C
! Note that this means we don't require a Fortran pre-precessor as the interface
! is a no-op when PAPI is not used.

module papi_interface

  use iso_c_binding

  implicit none

  interface

    subroutine papi_init() bind(C, name='papi_init')
    end subroutine

    subroutine papi_start() bind(C, name='papi_start')
    end subroutine

    subroutine papi_stop() bind(C, name='papi_stop')
    end subroutine

  end interface

end module papi_interface


