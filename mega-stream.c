
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
  const unsigned int elements = 100000000;
  const size_t size = sizeof(double)*elements;
  double *a = malloc(size);
  double *b = malloc(size);
  double *c = malloc(size);
  double *d = malloc(size);
  double *e = malloc(size);
  double *f = malloc(size);
  double *g = malloc(size);
  double *h = malloc(size);
  double *r = malloc(size);

  /* Initalise the data */
  #pragma omp parallel for
  for (unsigned int i = 0; i < elements; i++)
  {
    a[i] = 0.1;
    b[i] = 0.2;
    c[i] = 0.3;
    d[i] = 0.4;
    e[i] = 0.5;
    f[i] = 0.6;
    g[i] = 0.7;
    h[i] = 0.8;
    r[i] = 0.0;
  }


  for (int t = 0; t < 10; t++)
  {
    double tick = omp_get_wtime();
    /* Kernel */
    #pragma omp parallel for
    for (unsigned int i = 0; i < elements; i++)
    {
      r[i] = a[i]*b[i] + c[i]*d[i] + e[i]*f[i] + g[i]*h[i];
    }
    double tock = omp_get_wtime();
    printf("Iteration %d took %lf\n", t, tock-tick);
  }

  /* Check the results */
  double gold = 0.1*0.2 + 0.3*0.4 + 0.5*0.6 + 0.7*0.8;
  for (unsigned int i = 0; i < elements; i++)
  {
    if (r[i] != gold)
    {
      printf("Results incorrect\n");
      break;
    }
  }

  /* Free memory */
  free(a);
  free(b);
  free(c);
  free(d);
  free(e);
  free(f);
  free(g);
  free(h);
  free(r);

}
