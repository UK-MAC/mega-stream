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

#define VERSION "0.3"

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
#define IDX4(i,j,k,l,ni,nj,nk) ((i)+(ni)*IDX3((j),(k),(l),(nj),(nk)))
#define IDX5(i,j,k,l,m,ni,nj,nk,nl) ((i)+(ni)*IDX4((j),(k),(l),(m),(nj),(nk),(nl)))

/*
  Arrays are defined in terms of 3 sizes: inner, middle and outer.
  The large arrays are of size inner*middle*middle*middle*outer.
  The medium arrays are of size inner*middle*middle*outer.
  The small arrays are of size inner and are indexed with 1 index.

*/
#define OUTER   32 // 2^5
#define MIDDLE   8 // 2^3
#define INNER  128 // 2^7

/* Default alignment of 2 MB page boundaries */
#define ALIGNMENT 2*1024*1024

/* Tollerance with which to check final array values */
#define TOLR 1.0E-15

void parse_args(int argc, char *argv[]);

/* Default strides */
int Ni = INNER;
int Nj = MIDDLE;
int Nk = MIDDLE;
int Nl = MIDDLE;
int Nm = OUTER;

/* Number of iterations to run benchmark */
int ntimes = 100;

int main(int argc, char *argv[])
{

  printf("MEGA-STREAM! - v%s\n\n", VERSION);


  parse_args(argc, argv);

  printf("Small arrays:  %d elements\t\t(%.1lf KB)\n",
    Ni, Ni*sizeof(double)*1.0E-3);

  printf("Medium arrays: %d x %d x %d x %d elements\t(%.1lf MB)\n",
    Ni, Nj, Nj, Nm, Ni*Nj*Nj*Nm*sizeof(double)*1.0E-6);

  printf("Large arrays:  %d x %d x %d x %d x %d elements\t(%.1lf MB)\n",
    Ni, Nj, Nj, Nj, Nm, Ni*Nj*Nj*Nj*Nm*sizeof(double)*1.0E-6);

  const double footprint = (double)sizeof(double) * 1.0E-6 * (
    2.0*Ni*Nj*Nk*Nl*Nm +  /* r, q */
    Ni*Nj*Nk*Nm +         /* x */
    Ni*Nj*Nl*Nm +         /* y */
    Ni*Nk*Nl*Nm +         /* z */
    3.0*Ni +              /* a, b, c */
    Nj*Nk*Nl*Nm           /* sum */
    );
  printf("Memory footprint: %.1lf MB\n", footprint);

  /* Total memory moved - the arrays plus an extra sum as update is += */
  const double size = footprint + (double)sizeof(double) * L_size*M_size * 1.0E-6;

  printf("Running %d times\n", ntimes);

  printf("\n");

  double timings[ntimes];


  double *q = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);
  double *r = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);

  double *x = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nm);
  double *y = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nl*Nm);
  double *z = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nk*Nl*Nm);

  double *a = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);
  double *b = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);
  double *c = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);

  double *sum = aligned_alloc(ALIGNMENT, sizeof(double)*Nj*Nk*Nl*Nm);

  /* Initalise the data */
  #pragma omp parallel
  {
    /* q and r */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int l = 0; l < Nl; l++) {
        for (int k = 0; k < Nk; k++) {
          for (int j = 0; j < Nj; j++) {
            for (int i = 0; i < Ni; i++) {
              q[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] = 0.1;
              r[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] = 0.0;
            }
          }
        }
      }
    }

    /* x */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int k = 0; k < Nk; k++) {
        for (int j = 0; j < Nj; j++) {
          for (int i = 0; i < Ni; i++) {
            x[IDX4(i,j,k,m,Ni,Nj,Nk,Nm)] = 0.2;
          }
        }
      }
    }

    /* y */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int l = 0; l < Nl; l++) {
        for (int j = 0; j < Nj; j++) {
          for (int i = 0; i < Ni; i++) {
            y[IDX4(i,j,l,m,Ni,Nj,Nl,Nm)] = 0.3;
          }
        }
      }
    }

    /* z */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int l = 0; l < Nl; l++) {
        for (int k = 0; k < Nk; k++) {
          for (int i = 0; i < Ni; i++) {
            z[IDX4(i,k,l,m,Ni,Nk,Nl,Nm)] = 0.4;
          }
        }
      }
    }

    /* a, b, and c */
    #pragma omp for
    for (int i = 0; i < Ni; i++) {
      a[i] = 0.6;
      b[i] = 0.7;
      c[i] = 0.8;
    }

    /* sum */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int l = 0; l < Nl; l++) {
        for (int k = 0; k < Nk; k++) {
          for (int j = 0; j < Nj; j++) {
            sum[IDX4(j,k,l,m,Nj,Nk,Nl)] = 0.0;
          }
        }
      }
    }
  } /* End of parallel region */

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
