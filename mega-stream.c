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


  This aims to test the theory that streaming many large arrays causes memory
  bandwidth limits not to be reached, and latency becomes a dominating factor.
  We run a kernel with a similar form to the original triad, but with more than
  3 input arrays.

  The main kernel computes:
  r(i,j,k) = q(i,j,k) + a(i)*x(i,j) + b(i)*y(i,j) + c(i)*z(i,j)
  sum(j,k) = SUM(r(:,j,k))
*/

#define VERSION "0.2.1"

#include <float.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#define MAX(a,b) ((a) > (b)) ? (a) : (b)

#define IDX2(i,j,ni) ((i)+(ni)*(j))
#define IDX3(i,j,k,ni,nj) ((i)+(ni)*IDX2((j),(k),(nj)))

/*
  Arrays are defined in terms of 3 sizes
  The large arrays are of size SMALL*MEDIUM*LARGE and are indexed with 3 indicies.
  The medium arrays are of size SMALL*MEDIUM and are indexed with 2 indicies.
  The small arrays are of size SMALL and are indexed with 1 index.

  By default the large array has 2^27 elements, and the small array has 64 elements (2^6).
*/
#define LARGE  4096 // 2^12
#define MEDIUM  512 // 2^9
#define SMALL    64 // 2^6

/* Default alignment of 2 MB page boundaries */
#define ALIGNMENT 2*1024*1024

/* Tollerance with which to check final array values */
#define TOLR 1.0E-15

void parse_args(int argc, char *argv[]);

int L_size = LARGE;
int M_size = MEDIUM;
int S_size = SMALL;
int ntimes = 100;

int main(int argc, char *argv[])
{

  printf("MEGA-STREAM! - v%s\n\n", VERSION);


  parse_args(argc, argv);

  printf("Small arrays:  %d elements\t\t(%.1lf KB)\n",
    S_size, S_size*sizeof(double)*1.0E-3);

  printf("Medium arrays: %d x %d elements\t(%.1lf MB)\n",
    S_size, M_size, S_size*M_size*sizeof(double)*1.0E-6);

  printf("Large arrays:  %d x %d x %d elements\t(%.1lf MB)\n",
    S_size, M_size, L_size, S_size*M_size*L_size*sizeof(double)*1.0E-6);

  const double footprint = (double)sizeof(double) * 1.0E-6 * (
    2.0*L_size*M_size*S_size +   /* r, q */
    3.0*M_size*S_size +          /* x, y, z */
    3.0*S_size +                 /* a, b, c */
    L_size*M_size                /* sum */
    );
  printf("Memory footprint: %.1lf MB\n", footprint);

  /* Total memory moved - the arrays plus an extra sum as update is += */
  const double size = footprint + (double)sizeof(double) * L_size*M_size * 1.0E-6;

  printf("Running %d times\n", ntimes);

  printf("\n");

  double timings[ntimes];


  double *q = aligned_alloc(ALIGNMENT, sizeof(double)*L_size*M_size*S_size);
  double *r = aligned_alloc(ALIGNMENT, sizeof(double)*L_size*M_size*S_size);

  double *x = aligned_alloc(ALIGNMENT, sizeof(double)*M_size*S_size);
  double *y = aligned_alloc(ALIGNMENT, sizeof(double)*M_size*S_size);
  double *z = aligned_alloc(ALIGNMENT, sizeof(double)*M_size*S_size);

  double *a = aligned_alloc(ALIGNMENT, sizeof(double)*S_size);
  double *b = aligned_alloc(ALIGNMENT, sizeof(double)*S_size);
  double *c = aligned_alloc(ALIGNMENT, sizeof(double)*S_size);

  double *sum = aligned_alloc(ALIGNMENT, sizeof(double)*L_size*M_size);

  /* Initalise the data */
  #pragma omp parallel
  {
    #pragma omp for
    for (int k = 0; k < L_size; k++)
    {
      for (int j = 0; j < M_size; j++)
      {
        for (int i = 0; i < S_size; i++)
        {
          q[IDX3(i,j,k,S_size,M_size)] = 0.1;
          r[IDX3(i,j,k,S_size,M_size)] = 0.0;
        }
      }
    }

    #pragma omp for
    for (int j = 0; j < M_size; j++)
    {
      for (int i = 0; i < S_size; i++)
      {
        x[IDX2(i,j,S_size)] = 0.2;
        y[IDX2(i,j,S_size)] = 0.3;
        z[IDX2(i,j,S_size)] = 0.4;
      }
    }

    #pragma omp for
    for (int i = 0; i < S_size; i++)
    {
      a[i] = 0.6;
      b[i] = 0.7;
      c[i] = 0.8;
    }

    #pragma omp for
    for (int k = 0; k < L_size; k++)
    {
      for (int j = 0; j < M_size; j++)
      {
        sum[IDX2(j,k,M_size)] = 0.0;
      }
    }
  }

  /* Run the kernel multiple times */
  for (int t = 0; t < ntimes; t++)
  {
    double tick = omp_get_wtime();

    /**************************************************************************
     * Kernel
     *************************************************************************/
    #pragma omp parallel for
    for (int k = 0; k < L_size; k++)
    {
      for (int j = 0; j < M_size; j++)
      {
        double total = 0.0;
        #pragma omp simd reduction(+:total) aligned(r,q,a,b,c,x,y,z: ALIGNMENT)
        for (int i = 0; i < S_size; i++)
        {
          r[IDX3(i,j,k,S_size,M_size)] =
            q[IDX3(i,j,k,S_size,M_size)]
            + a[i] * x[IDX2(i,j,S_size)]
            + b[i] * y[IDX2(i,j,S_size)]
            + c[i] * z[IDX2(i,j,S_size)];

          total += r[IDX3(i,j,k,S_size,M_size)];
        }
        sum[IDX2(j,k,M_size)] += total;
      }
    }
    double tock = omp_get_wtime();
    timings[t] = tock-tick;

  }

  /* Check the results */
  const double gold = 0.1 + 0.2*0.6 + 0.3*0.7 + 0.4*0.8;
  const double gold_sum = gold*S_size*ntimes;

  /* Check the r array */
  for (int k = 0; k < L_size; k++)
    for (int j = 0; j < M_size; j++)
      for (int i = 0; i < S_size; i++)
      {
        if (fabs(r[IDX3(i,j,k,S_size,M_size)]-gold) > TOLR)
        {
          printf("Results incorrect - at (%d,%d,%d), %lf should be %lf\n",
            i,j,k, r[IDX3(i,j,k,S_size,M_size)], gold);
          goto sumcheck;
        }
      }

sumcheck:
  /* Check the reduction array */
  for (int i = 0; i < L_size*M_size; i++)
  {
    if (fabs(sum[i]-gold_sum) > TOLR)
    {
      printf("Reduction incorrect - at %d, %lf should be %lf\n",
        i, sum[i], gold_sum);
      break;
    }
  }

  /* Print timings */
  double min = DBL_MAX;
  double max = 0.0;
  double avg = 0.0;
  for (int t = 1; t < ntimes; t++)
  {
    min = MIN(min, timings[t]);
    max = MAX(max, timings[t]);
    avg += timings[t];
  }
  avg /= (double)(ntimes - 1);

  printf("Bandwidth MB/s  Min time    Max time    Avg time\n");
  printf("%12.1f %11.6f %11.6f %11.6f\n", size/min, min, max, avg);

  /* Free memory */
  free(q);
  free(r);
  free(x);
  free(y);
  free(z);
  free(a);
  free(b);
  free(c);
  free(sum);

  return EXIT_SUCCESS;

}

void parse_args(int argc, char *argv[])
{
  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "--large") == 0)
    {
      L_size = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--medium") == 0)
    {
      M_size = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--small") == 0)
    {
      S_size = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--swap") == 0)
    {
      int tmp = L_size;
      L_size = M_size;
      M_size = tmp;
    }
    else if (strcmp(argv[i], "--ntimes") == 0)
    {
      ntimes = atoi(argv[++i]);
      if (ntimes < 2)
      {
        fprintf(stderr, "ntimes must be 2 or greater\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(argv[i], "--help") == 0)
    {
      printf("Usage: %s [OPTION]\n", argv[0]);
      printf("\t --large n \tSet size of large dimension\n");
      printf("\t --medium n \tSet size of medium dimension\n");
      printf("\t --small n \tSet size of small dimension\n");
      printf("\t --swap\tSwap medium and large sizes over\n");
      printf("\t --ntimes n\tRun the benchmark n times\n");
      printf("\n");
      printf("\t Large  is %12d elements\n", LARGE);
      printf("\t Medium is %12d elements\n", MEDIUM);
      printf("\t Small  is %12d elements\n", SMALL);
      exit(EXIT_SUCCESS);
    }
    else
    {
      fprintf(stderr, "Unrecognised argument \"%s\"\n", argv[i]);
      exit(EXIT_FAILURE);
    }
  }
}
