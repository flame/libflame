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
#include "../src/libflame_interface.hh"
#define ENABLE_OMP_THREADS 1

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i,
      nrepeats, nb_algi, rval;

  double *Aref, *Ain, **Aout;
  int *piv_buff;

  FILE *fpin = fopen("dgetrf.in", "wb");
  FILE *fpout = fopen("dgetrf.out", "wb");

  nfirst = 80;
  nlast = 400;
  ninc = 80;
  nrepeats = 4;

  for ( n=nfirst; n<= nlast; n+=ninc ){

    Ain = (double *) malloc(n * n * sizeof(double));
    Aref = (double *) malloc(n * n * sizeof(double));
    piv_buff = (int *) malloc(n * sizeof(int));

    Aout = (double **) malloc(nrepeats * sizeof(double *));


    for(i = 0; i < nrepeats; i++)
    {
      Aout[i] = (double *) malloc(n * n * sizeof(double));
    }

    /* read input (generate random) */
    
    for(i = 0; i < n * n; i++)
    {
      Ain[i] = ( (double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    }

    for(i = 0; i < n * n; i++)
    {
      Aout[0][i] = Ain[i];
    }

    libflame::getrf(&n, &n, Aout[0], &n, piv_buff, &rval);

    {
       fwrite(Ain, sizeof(double), n * n, fpin);
       fwrite(Aout[0], sizeof(double), n * n, fpout);
       fwrite(piv_buff, sizeof(int), n, fpout);
    }

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
