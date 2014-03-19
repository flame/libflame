/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include <time.h>

#include "FLAME.h"
#include <omp.h>

#define N_VARIANTS 5

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_OPENMP_1TASK_FLA_PIP     1
#define FLA_ALG_OPENMP_1TASK_FLA_CIR     2
#define FLA_ALG_OPENMP_1TASK_REF_PIP     3
#define FLA_ALG_OPENMP_1TASK_REF_CIR     4


void time_Syrk_ln(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, n_input, m, n,
    nfirst, nlast, ninc,
    nb_alg,
    nrepeats,
    ntimings,
    variant,
    n_threads,
    n_thread_experiments,
    uplo,
    i, j,
    nvariants = N_VARIANTS;

  int n_threads_exp[64];

  char *colors = "brkgykkk";
  char *ticks = "o+*x-+*x";
  char matrix_config[3];
  char matrix_config_long[32];
  int  matrix_config_flag;

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff,
    d_n;

  FLA_Obj
    A, B, C, Cref;
  
  /* Initialize FLAME */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times */
  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &nrepeats );
  fprintf( stdout, "%c %d\n", '%', nrepeats );

  /* Blocking size used for blocked algorithms */
  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed */
  fprintf( stdout, "%c enter nfirst, nlast, ninc:", '%' );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  fprintf( stdout, "%c %d %d %d\n", '%', nfirst, nlast, ninc );

  /* Matrix dimensions */
  fprintf( stdout, "%c enter m (-1 means m=n, -m means swap m<=>n): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );
  
  fprintf( stdout, "%c enter number of thread experiments: ", '%' );
  scanf( "%d", &n_thread_experiments );
  fprintf( stdout, "%c %d\n", '%', n_thread_experiments );
	  
  fprintf( stdout, "%c enter number of threads for each experiment (separated by spaces): ", '%' );
  for( i = 0; i < n_thread_experiments; ++i )
    scanf( "%d", &n_threads_exp[i] );

  fprintf( stdout, "%c", '%' );
  for( i = 0; i < n_thread_experiments; ++i ) 
    fprintf( stdout, " %d", n_threads_exp[i] );


  /* Delete all existing data structures */
  fprintf( stdout, "\nclear all;\n\n" );


  if ( m_input == -1 ) 
  {
    matrix_config_flag = 0;
    strcpy(matrix_config,"sq");
    strcpy(matrix_config_long,"square matrix-matrix");
  }
  else if ( m_input < -1 )
  {
    matrix_config_flag = 1;
    strcpy(matrix_config,"op");
    strcpy(matrix_config_long,"outer panel-panel");
  }
  else
  {
    matrix_config_flag = 2;
    strcpy(matrix_config,"ip");
    strcpy(matrix_config_long,"inner panel-panel");
  }



  for ( n_input=nfirst; n_input<= nlast; n_input+=ninc ){
    
    if( matrix_config_flag == 0 ) 
    {
      m = n = n_input;
    }
    else if( matrix_config_flag == 1 )
    {
      m = n_input; 
      n = abs(m_input);
    }
    else
    {
      m = m_input; 
      n = n_input;
    }

    /* Allocate space for the matrices */

    /* Note: for some operations, there is no matrix B.  We are creating
       it anyway, for uniformity sake.
       Note: If your operation does something like B := A B than think of
       it as C := A C instead:  C is always the matrix that is the output. */
    FLA_Obj_create( FLA_DOUBLE, m, n, &A );
    FLA_Obj_create( FLA_DOUBLE, m, n, &B );
    FLA_Obj_create( FLA_DOUBLE, m, m, &C );
    FLA_Obj_create( FLA_DOUBLE, m, m, &Cref );

    /* Generate random matrices A, C */
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );
    FLA_Random_matrix( C );

    FLA_Copy_external( C, Cref );


    //FLA_Obj_show("C: ",C,"%f","\n");
    /* Time the reference implementation */
    time_Syrk_ln( 0, FLA_ALG_REFERENCE, nrepeats, n, nb_alg,
		A, B, C, Cref, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, (m_input > -1 ? n : m), gflops );
    fflush( stdout );

    for ( j = 0; j < n_thread_experiments; j++ ){

      n_threads = n_threads_exp[j];

      FLA_omp_set_num_threads( n_threads );
      FLA_omp_set_num_stages( n_threads );


      fprintf( stdout, "data_nth%d( %d, 1:9 ) = [ %d  ", n_threads, i, (m_input > -1 ? n : m) );
      fflush( stdout );

      time_Syrk_ln( 5, FLA_ALG_OPENMP_1TASK_FLA_PIP, nrepeats, n, nb_alg,
		  A, B, C, Cref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Syrk_ln( 5, FLA_ALG_OPENMP_1TASK_FLA_CIR, nrepeats, n, nb_alg,
		  A, B, C, Cref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Syrk_ln( 5, FLA_ALG_OPENMP_1TASK_REF_PIP, nrepeats, n, nb_alg,
		  A, B, C, Cref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Syrk_ln( 5, FLA_ALG_OPENMP_1TASK_REF_CIR, nrepeats, n, nb_alg,
		  A, B, C, Cref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &Cref );
    fprintf( stdout, "\n" );


  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  fprintf( stdout, "figure;\n" );

  /* Plot the performance of the reference implementation */
  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on;\n" );

  /* Plot the data for the other numbers of threads */
  for ( i = 0; i < n_thread_experiments; i++ ){
    fprintf( stdout, "plot( data_nth%d( :,1 ), data_nth%d( :, 2 ), '%c-%c' ); \n", 
                      n_threads_exp[ i ], n_threads_exp[ i ], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_nth%d( :,1 ), data_nth%d( :, 4 ), '%c--%c' ); \n", 
                      n_threads_exp[ i ], n_threads_exp[ i ], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_nth%d( :,1 ), data_nth%d( :, 6 ), '%c-.%c' ); \n", 
                      n_threads_exp[ i ], n_threads_exp[ i ], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_nth%d( :,1 ), data_nth%d( :, 8 ), '%c:%c' ); \n", 
                      n_threads_exp[ i ], n_threads_exp[ i ], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( 'Reference', ... \n" );

  for ( i = 0; i < n_thread_experiments-1; i++ )
    printf( "'fp\\_nth%d', 'fc\\_nth%d', 'rp\\_nth%d', 'rc\\_nth%d', ... \n", 
             n_threads_exp[ i ], 
             n_threads_exp[ i ], 
             n_threads_exp[ i ], 
             n_threads_exp[ i ] );

  printf( "'fp\\_nth%d', 'fc\\_nth%d', 'rp\\_nth%d', 'rc\\_nth%d', 2 );\n", 
           n_threads_exp[ n_thread_experiments-1 ], 
           n_threads_exp[ n_thread_experiments-1 ], 
           n_threads_exp[ n_thread_experiments-1 ], 
           n_threads_exp[ n_thread_experiments-1 ] );

  fprintf( stdout, "xlabel( 'matrix dimension %c' );\n", (m_input > -1 ? 'n' : 'm') );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n");
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", nlast, n_threads_exp[n_thread_experiments-1] * max_gflops );
  fprintf( stdout, "title( 'OpenFLAME syrk\\_ln\\_var5 performance (one task, %s)' );\n", matrix_config_long );
  fprintf( stdout, "print -depsc omp1t_var5_%s.eps\n", matrix_config);
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

