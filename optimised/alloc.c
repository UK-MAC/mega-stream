/*

  Copyright 2016 Tom Deakin, University of Bristol

  This file is part of mega-stream.

  mega-stream is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  mega-stream is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with mega-stream.  If not, see <http://www.gnu.org/licenses/>.


  This aims to investigate the limiting factor for a simple kernel, in particular
  where bandwidth limits not to be reached, and latency becomes a dominating factor.

*/

#include <stdlib.h>
#include <stdio.h>

#define ALIGNMENT 2*1024*1024

double * alloc(int *len) {
  double * p = (double *)aligned_alloc(ALIGNMENT, sizeof(double)*(*len));
  printf("Allocated %p size %zu. len was %d\n", p, sizeof(double)*(*len), (*len));
  return p;
}

void alloc_free(double ** p)  {
  printf("About to free %p\n", *p);
  free(*p);
}

