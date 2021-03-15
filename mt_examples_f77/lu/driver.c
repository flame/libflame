/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "omp.h"

#define ENABLE_OMP_THREADS 1

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, j,
      nrepeats, nb_algi, rval;

  double *Aref, *Ain, **Aout;
  int **piv_buff, *piv_buff_ref;

  double max_diff, avg_max_diff;

  FILE *fpin = fopen("dgetrf.in", "rb");
  FILE *fpout = fopen("dgetrf.out", "rb");

  nfirst = 80;
  nlast = 400;
  ninc = 80;
  nrepeats = 4;

  for(n = nfirst; n <= nlast; n += ninc)
  {
    Ain = (double *) malloc(n * n * sizeof(double));
    Aref = (double *) malloc(n * n * sizeof(double));
    piv_buff_ref = (int *) malloc(n * sizeof(int));

    piv_buff = (int **) malloc(nrepeats * sizeof(int *));
    Aout = (double **) malloc(nrepeats * sizeof(double *));


    for(i = 0; i < nrepeats; i++)
    {
      piv_buff[i] = (int *) malloc(n * sizeof(int));
      Aout[i] = (double *) malloc(n * n * sizeof(double));
    }
    Aout[0][0] = 0;
    /* read input and ref output (generate random) */
    fread(Ain, sizeof(double), n * n, fpin);
    fread(Aref, sizeof(double), n * n, fpout);
    fread(piv_buff_ref, sizeof(int), n, fpout);
    
    /* Initialize the output with input and send it to lapack function */
    for(i = 0; i < nrepeats; i++)
    {
      for(j = 0; j < n * n; j++)
      {
        Aout[i][j] = Ain[j];
      }
    }
#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for(i = 0; i < nrepeats; i++)
    {
      dgetrf_(&n, &n, Aout[i], &n, piv_buff[i], &rval);
    }

    avg_max_diff = 0;
    for(i = 0; i < nrepeats; i++)
    {
      /* Max. Average difference */
      max_diff = 0;
      for(j = 0; j < n * n; j++)
      {
        max_diff = fmax(abs(Aout[i][j] - Aref[j]), max_diff);
      }
      avg_max_diff += (max_diff / nrepeats);
    }

    printf("Matrix Size:%dx%d\nMax Difference:%lf\n\n", n, n, avg_max_diff);

    free(Ain);
    free(Aref);
    free(piv_buff);

    for(i = 0; i < nrepeats; i++)
    {
      free(Aout[i]);
    }

    free(Aout);
  }

  fclose(fpin);
  fclose(fpout);
}
