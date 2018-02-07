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

#define VERSION "2.0"

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
#define OUTER   64 // 3^6
#define MIDDLE  16 // 2^4
#define INNER  128 // 2^7

/* Default alignment of 2 MB page boundaries */
#define ALIGNMENT 2*1024*1024

/* Tollerance with which to check final array values */
#define TOLR 1.0E-15

/* Starting values */
#define R_START 0.0
#define Q_START 0.01
#define X_START 0.02
#define Y_START 0.03
#define Z_START 0.04
#define A_START 0.06
#define B_START 0.07
#define C_START 0.08


void kernel(
  const int Ni, const int Nj, const int Nk, const int Nl, const int Nm,
  double * restrict r, const double * restrict q,
  double * restrict x, double * restrict y, double * restrict z,
  const double * restrict a, const double * restrict b, const double * restrict c,
  double * restrict sum
);
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

  printf("Medium arrays: %d x %d x %d x %d elements\t(%.1lf MB)\n"
         "               %d x %d x %d x %d elements\t(%.1lf MB)\n"
         "               %d x %d x %d x %d elements\t(%.1lf MB)\n",
    Ni, Nj, Nk, Nm, Ni*Nj*Nk*Nm*sizeof(double)*1.0E-6,
    Ni, Nj, Nl, Nm, Ni*Nj*Nl*Nm*sizeof(double)*1.0E-6,
    Ni, Nk, Nl, Nm, Ni*Nk*Nl*Nm*sizeof(double)*1.0E-6);

  printf("Large arrays:  %d x %d x %d x %d x %d elements\t(%.1lf MB)\n",
    Ni, Nj, Nk, Nl, Nm, Ni*Nj*Nk*Nl*Nm*sizeof(double)*1.0E-6);

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
  const double moved = (double)sizeof(double) * 1.0E-6 * (
    Ni*Nj*Nk*Nl*Nm  + /* read q */
    Ni*Nj*Nk*Nl*Nm  + /* write r */
    Ni + Ni + Ni    + /* read a, b and c */
    2.0*Ni*Nj*Nk*Nm + /* read and write x */
    2.0*Ni*Nj*Nl*Nm + /* read and write y */
    2.0*Ni*Nk*Nl*Nm + /* read and write z */
    2.0*Nj*Nk*Nl*Nm   /* read and write sum */
  );

  printf("Running %d times\n", ntimes);

  printf("\n");

  double timings[ntimes];

#ifdef __APPLE__
  double *q;
  posix_memalign((void **)&q, ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);
  double *r;
  posix_memalign((void **)&r, ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);

  double *x;
  posix_memalign((void **)&x, ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nm);
  double *y;
  posix_memalign((void **)&y, ALIGNMENT, sizeof(double)*Ni*Nj*Nl*Nm);
  double *z;
  posix_memalign((void **)&z, ALIGNMENT, sizeof(double)*Ni*Nk*Nl*Nm);

  double *a;
  posix_memalign((void **)&a, ALIGNMENT, sizeof(double)*Ni);
  double *b;
  posix_memalign((void **)&b, ALIGNMENT, sizeof(double)*Ni);
  double *c;
  posix_memalign((void **)&c, ALIGNMENT, sizeof(double)*Ni);

  double *sum;
  posix_memalign((void **)&sum, ALIGNMENT, sizeof(double)*Nj*Nk*Nl*Nm);
#else
  double *q = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);
  double *r = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nl*Nm);

  double *x = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nk*Nm);
  double *y = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nj*Nl*Nm);
  double *z = aligned_alloc(ALIGNMENT, sizeof(double)*Ni*Nk*Nl*Nm);

  double *a = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);
  double *b = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);
  double *c = aligned_alloc(ALIGNMENT, sizeof(double)*Ni);

  double *sum = aligned_alloc(ALIGNMENT, sizeof(double)*Nj*Nk*Nl*Nm);
#endif

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
              q[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] = Q_START;
              r[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] = R_START;
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
            x[IDX4(i,j,k,m,Ni,Nj,Nk)] = X_START;
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
            y[IDX4(i,j,l,m,Ni,Nj,Nl)] = Y_START;
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
            z[IDX4(i,k,l,m,Ni,Nk,Nl)] = Z_START;
          }
        }
      }
    }

    /* a, b, and c */
    #pragma omp for
    for (int i = 0; i < Ni; i++) {
      a[i] = A_START;
      b[i] = B_START;
      c[i] = C_START;
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

  double begin = omp_get_wtime();
  /* Run the kernel multiple times */
  for (int t = 0; t < ntimes; t++) {
    double tick = omp_get_wtime();

    kernel(Ni,Nj,Nk,Nl,Nm,r,q,x,y,z,a,b,c,sum);

    /* Swap the pointers */
    double *tmp = q; q = r; r = tmp;

    double tock = omp_get_wtime();
    timings[t] = tock-tick;

  }

  double end = omp_get_wtime();

  /* Check the results - total of the sum array */
  double total = 0.0;
  for (int i = 0; i < Nj*Nk*Nl*Nm; i++)
    total += sum[i];
  printf("Sum total: %lf\n", total);

  /* Print timings */
  double min = DBL_MAX;
  double max = 0.0;
  double avg = 0.0;
  for (int t = 1; t < ntimes; t++) {
    min = MIN(min, timings[t]);
    max = MAX(max, timings[t]);
    avg += timings[t];
  }
  avg /= (double)(ntimes - 1);

  printf("\n");
  printf("Bandwidth MB/s  Min time    Max time    Avg time\n");
  printf("%12.1f %11.6f %11.6f %11.6f\n", moved/min, min, max, avg);
  printf("Total time: %11.6f\n", end-begin);

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

/**************************************************************************
 * Kernel
 *************************************************************************/
void kernel(
  const int Ni, const int Nj, const int Nk, const int Nl, const int Nm,
  double * restrict r,
  const double * restrict q,
  double * restrict x,
  double * restrict y,
  double * restrict z,
  const double * restrict a,
  const double * restrict b,
  const double * restrict c,
  double * restrict sum
  )
{
  #pragma omp parallel for
  for (int m = 0; m < Nm; m++) {
    for (int l = 0; l < Nl; l++) {
      for (int k = 0; k < Nk; k++) {
        for (int j = 0; j < Nj; j++) {
          double total = 0.0;
          #pragma omp simd reduction(+:total)
          for (int i = 0; i < Ni; i++) {
            /* Set r */
            double tmp_r =
              q[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] +
              a[i] * x[IDX4(i,j,k,m,Ni,Nj,Nk)] +
              b[i] * y[IDX4(i,j,l,m,Ni,Nj,Nl)] +
              c[i] * z[IDX4(i,k,l,m,Ni,Nk,Nl)];

            /* Update x, y and z */
            x[IDX4(i,j,k,m,Ni,Nj,Nk)] = 0.2*tmp_r - x[IDX4(i,j,k,m,Ni,Nj,Nk)];
            y[IDX4(i,j,l,m,Ni,Nj,Nl)] = 0.2*tmp_r - y[IDX4(i,j,l,m,Ni,Nj,Nl)];
            z[IDX4(i,k,l,m,Ni,Nk,Nl)] = 0.2*tmp_r - z[IDX4(i,k,l,m,Ni,Nk,Nl)];

            /* Reduce over Ni */
            total += tmp_r;

            /* Save r */
            r[IDX5(i,j,k,l,m,Ni,Nj,Nk,Nl)] = tmp_r;

          } /* Ni */

          sum[IDX4(j,k,l,m,Nj,Nk,Nl)] += total;

        } /* Nj */
      } /* Nk */
    } /* Nl */
  } /* Nm */
}

void parse_args(int argc, char *argv[])
{
  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "--outer") == 0)
    {
      Nm = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--inner") == 0)
    {
      Ni = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--middle") == 0)
    {
      int num = atoi(argv[++i]);
      Nj = num;
      Nk = num;
      Nl = num;
    }
    else if (strcmp(argv[i], "--Nj") == 0)
    {
      Nj = atoi(argv[++i]);
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
      printf("\t --outer  n \tSet size of outer dimension\n");
      printf("\t --inner  n \tSet size of middle dimensions\n");
      printf("\t --middle n \tSet size of inner dimension\n");
      printf("\t --Nj     n \tSet size of the j dimension\n");
      printf("\t --ntimes n\tRun the benchmark n times\n");
      printf("\n");
      printf("\t Outer   is %12d elements\n", OUTER);
      printf("\t Middle are %12d elements\n", MIDDLE);
      printf("\t Inner   is %12d elements\n", INNER);
      exit(EXIT_SUCCESS);
    }
    else
    {
      fprintf(stderr, "Unrecognised argument \"%s\"\n", argv[i]);
      exit(EXIT_FAILURE);
    }
  }
}
