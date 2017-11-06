#
# Copyright 2016 Tom Deakin, University of Bristol
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
FFLAGS_ = -O3

# Intel settings
FTN_INTEL = mpiifort
FFLAGS_INTEL = -O3 -xHost -align array256byte -qopt-report=5
OMP_INTEL = -qopenmp

# GNU settings
FTN_GNU = mpif90
FFLAGS_GNU = -O3 -ffree-line-length-none
OMP_GNU = -fopenmp

# Cray settings
FTN_CRAY = ftn
FFLAGS_CRAY = -O3
OMP_CRAY =

ifdef DEBUG
  FFLAGS_GNU = -ffree-line-length-none -g -fcheck=all
  FFLAGS_INTEL = -g -debug all -traceback
endif

# Choose compiler
FTN=$(FTN_$(COMPILER))
FFLAGS=$(FFLAGS_$(COMPILER))
OMP=$(OMP_$(COMPILER))

default: mega-sweep

mega-sweep: comms.f90 kernel.f90 kernely.f90 mega-sweep.f90
	$(FTN) $(FFLAGS) $(OPTIONS) $(OMP) $^ -o $@

.PHONY: clean
clean:
	rm -f mega-sweep comms.mod

