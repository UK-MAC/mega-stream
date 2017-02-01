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
#include <assert.h>

#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#define MAX(a,b) ((a) > (b)) ? (a) : (b)

#define IDX2(i,j,ni) ((i)+(ni)*(j))
#define IDX3(i,j,k,ni,nj) ((i)+(ni)*IDX2((j),(k),(nj)))
#define IDX4(i,j,k,l,ni,nj,nk) ((i)+(ni)*IDX3((j),(k),(l),(nj),(nk)))
#define IDX5(i,j,k,l,m,ni,nj,nk,nl) ((i)+(ni)*IDX4((j),(k),(l),(m),(nj),(nk),(nl)))
#define IDX6(i,j,k,l,m,n,ni,nj,nk,nl,nm) ((i)+(ni)*IDX5((j),(k),(l),(m),(n),(nj),(nk),(nl),(nm)))

/*
  Arrays are defined in terms of 3 sizes: inner, middle and outer.
  The large arrays are of size inner*middle*middle*middle*outer.
  The medium arrays are of size inner*middle*middle*outer.
  The small arrays are of size inner and are indexed with 1 index.

*/
#define OUTER   64 // 3^6
#define MIDDLE  16 // 2^4
#define INNER  128 // 2^7

/* Vector length is machine-specific and should be overridden by makefile */
#ifndef VLEN
#define VLEN 8
#endif
static_assert((VLEN > 0) && ((VLEN & (VLEN-1)) == 0), "VLEN must be a power of 2.");

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

#ifdef __APPLE__
void* aligned_alloc(size_t alignment, size_t size)
{
    void* mem;
    posix_memalign(&mem, alignment, size);
    return mem;
}
#endif

void kernel(
  const int Ng,
  const int Ni, const int Nj, const int Nk, const int Nl, const int Nm,
  double (* restrict r)[Ng][Nl][Nk][Nj][VLEN],
  const double (* restrict q)[Ng][Nl][Nk][Nj][VLEN],
  double (* restrict x)[Ng][Nk][Nj][VLEN],
  double (* restrict y)[Ng][Nl][Nj][VLEN],
  double (* restrict z)[Ng][Nl][Nk][VLEN],
  const double (* restrict a)[VLEN],
  const double (* restrict b)[VLEN],
  const double (* restrict c)[VLEN],
  double (* restrict sum)[Nl][Nk][Nj]
);
void parse_args(int argc, char *argv[]);

/* Default strides */
int Ni = INNER;
int Nj = MIDDLE;
int Nk = MIDDLE;
int Nl = MIDDLE;
int Nm = OUTER;
int Ng;

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

  /* Total memory moved (in the best case) - the arrays plus an extra sum as update is += */
  const double moved = (double)sizeof(double) * 1.0E-6 * (
    Ni*Nj*Nk*Nl*Nm  + /* read q */
    Ni*Nj*Nk*Nl*Nm  + /* write r */
    Ni + Ni + Ni    + /* read a, b and c */
    2.0*Ni*Nj*Nk*Nm + /* read and write x */
    2.0*Ni*Nj*Nl*Nm + /* read and write y */
    2.0*Ni*Nk*Nl*Nm + /* read and write z */
    2.0*Nj*Nk*Nl*Nm   /* read and write sum */
  );

  /* Split inner-most dimension into VLEN-sized chunks */
  Ng = ((Ni + (VLEN-1)) & ~(VLEN-1)) / VLEN;
  printf("Inner dimension split into %d chunks of size %d\n", Ng, VLEN);

  printf("Running %d times\n", ntimes);

  printf("\n");

  double timings[ntimes];

  double *q = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Nj*Nk*Nl*Nm*Ng);
  double *r = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Nj*Nk*Nl*Nm*Ng);

  double *x = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Nj*Nk*Nm*Ng);
  double *y = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Nj*Nl*Nm*Ng);
  double *z = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Nk*Nl*Nm*Ng);

  double *a = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Ng);
  double *b = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Ng);
  double *c = aligned_alloc(ALIGNMENT, sizeof(double)*VLEN*Ng);

  double *sum = aligned_alloc(ALIGNMENT, sizeof(double)*Nj*Nk*Nl*Nm);

  /* Initalise the data */
  #pragma omp parallel
  {
    /* q and r */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int g = 0; g < Ng; g++) {
        for (int l = 0; l < Nl; l++) {
          for (int k = 0; k < Nk; k++) {
            for (int j = 0; j < Nj; j++) {
              #pragma omp simd
              for (int v = 0; v < VLEN; v++) {
                q[IDX6(v,j,k,l,g,m,VLEN,Nj,Nk,Nl,Ng)] = Q_START;
                r[IDX6(v,j,k,l,g,m,VLEN,Nj,Nk,Nl,Ng)] = R_START;
              }
            }
          }
        }
      }
    }

    /* x */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int g = 0; g < Ng; g++) {
        for (int k = 0; k < Nk; k++) {
          for (int j = 0; j < Nj; j++) {
            #pragma omp simd
            for (int v = 0; v < VLEN; v++) {
              x[IDX5(v,j,k,g,m,VLEN,Nj,Nk,Ng)] = X_START;
            }
          }
        }
      }
    }

    /* y */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int g = 0; g < Ng; g++) {
        for (int l = 0; l < Nl; l++) {
          for (int j = 0; j < Nj; j++) {
            #pragma omp simd
            for (int v = 0; v < VLEN; v++) {
              y[IDX5(v,j,l,g,m,VLEN,Nj,Nl,Ng)] = Y_START;
            }
          }
        }
      }
    }

    /* z */
    #pragma omp for
    for (int m = 0; m < Nm; m++) {
      for (int g = 0; g < Ng; g++) {
        for (int l = 0; l < Nl; l++) {
          for (int k = 0; k < Nk; k++) {
            #pragma omp simd
            for (int v = 0; v < VLEN; v++) {
              z[IDX5(v,k,l,g,m,VLEN,Nk,Nl,Ng)] = Z_START;
            }
          }
        }
      }
    }

    /* a, b, and c */
    #pragma omp for
    for (int g = 0; g < Ng; g++) {
      #pragma omp simd
      for (int v = 0; v < VLEN; v++) {
        a[IDX2(v,g,VLEN)] = A_START;
        b[IDX2(v,g,VLEN)] = B_START;
        c[IDX2(v,g,VLEN)] = C_START;
      }
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
  for (int t = 0; t < ntimes; t++) {
    double tick = omp_get_wtime();

    kernel(Ng,Ni,Nj,Nk,Nl,Nm,
        (double(*)[Ng][Nl][Nk][Nj][VLEN]) r,
        (const double(*)[Ng][Nl][Nk][Nj][VLEN]) q,
        (double(*)[Ng][Nk][Nj][VLEN]) x,
        (double(*)[Ng][Nl][Nj][VLEN]) y,
        (double(*)[Ng][Nl][Nk][VLEN]) z,
        (const double(*)[VLEN]) a,
        (const double(*)[VLEN]) b,
        (const double(*)[VLEN]) c,
        (double(*)[Nl][Nk][Nj]) sum);

    /* Swap the pointers */
    double *tmp = q; q = r; r = tmp;

    double tock = omp_get_wtime();
    timings[t] = tock-tick;

  }

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
#include <immintrin.h>
void kernel(
  const int Ng,
  const int Ni, const int Nj, const int Nk, const int Nl, const int Nm,
  double (* restrict r)[Ng][Nl][Nk][Nj][VLEN],
  const double (* restrict q)[Ng][Nl][Nk][Nj][VLEN],
  double (* restrict x)[Ng][Nk][Nj][VLEN],
  double (* restrict y)[Ng][Nl][Nj][VLEN],
  double (* restrict z)[Ng][Nl][Nk][VLEN],
  const double (* restrict a)[VLEN],
  const double (* restrict b)[VLEN],
  const double (* restrict c)[VLEN],
  double (* restrict sum)[Nl][Nk][Nj]
  )
{
  #pragma omp parallel for
  for (int m = 0; m < Nm; m++) {
    for (int g = 0; g < Ng; g++) {
      for (int l = 0; l < Nl; l++) {
        for (int k = 0; k < Nk; k++) {
          for (int j = 0; j < Nj; j++) {
            double total = 0.0;
            _mm_prefetch((const char*) (&q[m][g][l][k][j][0] + 32*VLEN), _MM_HINT_T1);
            #pragma vector nontemporal(r)
            #pragma omp simd reduction(+:total) aligned(a,b,c,x,y,z,r,q:64)
            for (int v = 0; v < VLEN; v++) {
              /* Set r */
              r[m][g][l][k][j][v] =
                q[m][g][l][k][j][v] +
                a[g][v] * x[m][g][k][j][v] +
                b[g][v] * y[m][g][l][j][v] +
                c[g][v] * z[m][g][l][k][v];

              /* Update x, y and z */
              x[m][g][k][j][v] = 0.2*r[m][g][l][k][j][v] - x[m][g][k][j][v];
              y[m][g][l][j][v] = 0.2*r[m][g][l][k][j][v] - y[m][g][l][j][v];
              z[m][g][l][k][v] = 0.2*r[m][g][l][k][j][v] - z[m][g][l][k][v];

              /* Reduce over Ni */
              total += r[m][g][l][k][j][v];

            } /* VLEN */

            sum[m][l][k][j] += total;

          } /* Nj */
        } /* Nk */
      } /* Nl */
    } /* Ng */
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
