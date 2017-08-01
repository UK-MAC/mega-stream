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


FTN = mpiifort
FFLAGS = -O3 -xHost -align array256byte -qopt-report=5

OMP = -qopenmp

default: mega-sweep

mega-sweep: mega-sweep.f90 comms.f90
	$(FTN) $(FFLAGS) $(OPTIONS) $(OMP) $^ -o $@

.PHONY: clean
clean:
	rm -f mega-sweep

