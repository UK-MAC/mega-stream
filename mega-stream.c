
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
*/

#define VERSION "0.1"

#include <float.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>

#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#define MAX(a,b) ((a) > (b)) ? (a) : (b)

#define LARGE   134217728// 2^27
#define MEDIUM    8388608 // 2^23
#define SMALL         128

void parse_args(int argc, char *argv[]);

unsigned int L_size = LARGE;
unsigned int M_size = MEDIUM;
unsigned int S_size = SMALL;

int main(int argc, char *argv[])
{

  printf("MEGA-STREAM! - v%s\n", VERSION);

  parse_args(argc, argv);

  const int ntimes = 100;
  double timings[ntimes];

  const double size = 8.0 * (2.0*L_size + 3.0*M_size + 3.0*S_size) * 1.0E-6;

  /* The following assumes sizes are powers of 2 */
  const unsigned int S_mask = S_size - 1;
  const unsigned int M_mask = M_size - 1;

  /* Vectorization assumes that sizes are multiples of 8 */
  assert(S_size % 8 == 0);
  assert(M_size % 8 == 0);
  assert(L_size % 8 == 0);

  /* Align arrays to 2 MB page boundaries */
  const unsigned int alignment = 2*1024*1024;
  double *restrict const q = _mm_malloc(sizeof(double)*L_size, alignment);
  double *restrict const r = _mm_malloc(sizeof(double)*L_size, alignment);

  double *restrict const x = _mm_malloc(sizeof(double)*M_size, alignment);
  double *restrict const y = _mm_malloc(sizeof(double)*M_size, alignment);
  double *restrict const z = _mm_malloc(sizeof(double)*M_size, alignment);

  double *restrict const a = _mm_malloc(sizeof(double)*S_size, alignment);
  double *restrict const b = _mm_malloc(sizeof(double)*S_size, alignment);
  double *restrict const c = _mm_malloc(sizeof(double)*S_size, alignment);

  /* Initalise the data */
  #pragma omp parallel
  {
    #pragma omp for
    for (int i = 0; i < L_size; i++)
    {
      q[i] = 0.1;
      r[i] = 0.0;
    }

    #pragma omp for
    for (int i = 0; i < M_size; i++)
    {
      x[i] = 0.2;
      y[i] = 0.3;
      z[i] = 0.4;
    }

    #pragma omp for
    for (int i = 0; i < S_size; i++)
    {
      a[i] = 0.6;
      b[i] = 0.7;
      c[i] = 0.8;
    }

  }

  /* Run the kernel multiple times */
  for (int t = 0; t < ntimes; t++)
  {
    double tick = omp_get_wtime();
    /* Kernel */
    #pragma omp parallel for
    for (int i = 0; i < L_size; i += 8)
    {
      int s = i&S_mask;
      int m = i&M_mask;
      #pragma vector aligned
      #pragma vector nontemporal
      #pragma omp simd
      for (int v = 0; v < 8; ++v)
      {
         r[i+v] = q[i+v] + a[s+v]*x[m+v] + b[s+v]*y[m+v] + c[s+v]*z[m+v];
      }
    }
    double tock = omp_get_wtime();
    timings[t] = tock-tick;
  }

  /* Check the results */
  double gold = 0.1 + 0.2*0.6 + 0.3*0.7 + 0.4*0.8;
  for (int i = 0; i < L_size; i++)
  {
    if (r[i] != gold)
    {
      printf("Results incorrect\n");
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
  _mm_free(q);
  _mm_free(r);
  _mm_free(x);
  _mm_free(y);
  _mm_free(z);
  _mm_free(a);
  _mm_free(b);
  _mm_free(c);

  return EXIT_SUCCESS;

}

void parse_args(int argc, char *argv[])
{
  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "--large") == 0)
    {
      /* All arrays are large */
      printf("Setting: Large\n");
      M_size = LARGE;
      S_size = LARGE;
    }
    else if (strcmp(argv[i], "--medium") == 0)
    {
      /* Large arrays with only medium arrays */
      printf("Setting: Medium\n");
      S_size = MEDIUM;
    }
    else if (strcmp(argv[i], "--small") == 0)
    {
      /* Large arrays with only small arrays */
      printf("Setting: Small\n");
      M_size = SMALL;
    }
    else if (strcmp(argv[i], "--custom") == 0)
    {
      unsigned int size = atoi(argv[++i]);
      printf("Setting: Custom - %d\n", size);
      M_size = size;
      S_size = size;
    }
    else if (strcmp(argv[i], "--help") == 0)
    {
      printf("Usage: %s [OPTION]\n", argv[0]);
      printf("\t --large\tMake all arrays large in size\n");
      printf("\t --medium\t2 large arrays, and 6 medium arrays\n");
      printf("\t --small\t2 large arrays, and 6 small arrays\n");
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
