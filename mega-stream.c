
/*
  This aims to test the theory that streaming many large arrays causes memory
  bandwidth limits not to be reached, and latency becomes a dominating factor.
  We run a kernel with a similar form to the original triad, but with more than
  3 input arrays.
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{

  const int ntimes = 10;

  const unsigned int L_size = 200000000;
  const unsigned int M_size =  10000000;
  const unsigned int S_size =       100;

  const double size = 8.0 * (2.0*L_size + 3.0*M_size + 3.0*S_size) * 1.0E-6;

  double *q = malloc(sizeof(double)*L_size);
  double *r = malloc(sizeof(double)*L_size);

  double *x = malloc(sizeof(double)*M_size);
  double *y = malloc(sizeof(double)*M_size);
  double *z = malloc(sizeof(double)*M_size);

  double *a = malloc(sizeof(double)*S_size);
  double *b = malloc(sizeof(double)*S_size);
  double *c = malloc(sizeof(double)*S_size);

  /* Initalise the data */
  #pragma omp parallel
  {
    #pragma omp for
    for (unsigned int i = 0; i < L_size; i++)
    {
      q[i] = 0.1;
      r[i] = 0.0;
    }

    #pragma omp for
    for (unsigned int i = 0; i < M_size; i++)
    {
      x[i] = 0.2;
      y[i] = 0.3;
      z[i] = 0.4;
    }

    #pragma omp for
    for (unsigned int i = 0; i < S_size; i++)
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
    for (unsigned int i = 0; i < L_size; i++)
    {
      r[i] = q[i] + a[i%S_size]*x[i%M_size] + b[i%S_size]*y[i%M_size] + c[i%S_size]*z[i%M_size];
    }
    double tock = omp_get_wtime();
    printf("Iteration %d took %lf\n", t, tock-tick);
    printf("Bandwidth: %lf MB/s\n", size/(tock-tick));

  }

  /* Check the results */
  double gold = 0.1 + 0.2*0.6 + 0.3*0.7 + 0.4*0.8;
  for (unsigned int i = 0; i < L_size; i++)
  {
    if (r[i] != gold)
    {
      printf("Results incorrect\n");
      break;
    }
  }

  /* Free memory */
  free(q);
  free(r);
  free(x);
  free(y);
  free(z);
  free(a);
  free(b);
  free(c);

}
