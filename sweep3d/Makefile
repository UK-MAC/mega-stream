#
# Copyright 2017 Tom Deakin, University of Bristol
#
# This file is part of mega-stream.
#
# mega-stream is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mega-stream is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mega-stream.  If not, see <http://www.gnu.org/licenses/>.
#

# Default compiler and flags
FTN_ = mpif90
CC_ = cc
FFLAGS_ = -O3

# Intel settings
FTN_INTEL = mpiifort
CC_INTEL = icc
FFLAGS_INTEL = -O3 -xHost -align array256byte -qopt-report=5
OMP_INTEL = -qopenmp

# GNU settings
FTN_GNU = mpif90
CC_GNU = gcc
FFLAGS_GNU = -O3 -ffree-line-length-none
OMP_GNU = -fopenmp

# Cray settings
FTN_CRAY = ftn
CC_CRAY = cc
FFLAGS_CRAY = -O3
OMP_CRAY =

ifdef DEBUG
  FFLAGS_GNU = -ffree-line-length-none -g -fcheck=all
  FFLAGS_INTEL = -g -debug all -traceback
endif

# Choose compiler
FTN=$(FTN_$(COMPILER))
CC=$(CC_$(COMPILER))
FFLAGS=$(FFLAGS_$(COMPILER))
OMP=$(OMP_$(COMPILER))

# PAPI library
ifdef PAPI
  PAPI_LIB = -DPAPI -lpapi
endif

default: mega-sweep3d

papi.o papi_interface.mod: papi.c
	$(CC) papi.c -c $(OMP) $(PAPI_LIB)
	# Ensure the .mod file is produced before we try to compile
	$(FTN) papi_interface.f90 -c

mega-sweep3d: comms3d.f90 kernel3d.f90 mega-sweep3d.f90 papi_interface.f90 papi.o
	$(FTN) $(FFLAGS) $(OPTIONS) $(OMP) $^ -o $@ $(PAPI_LIB)

.PHONY: clean
clean:
	rm -f mega-sweep comms.mod

