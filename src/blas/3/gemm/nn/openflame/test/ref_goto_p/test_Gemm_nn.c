/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include <time.h>

#include "FLAME.h"

#define N_VARIANTS 5

#define FLA_ALG_REFERENCE 0

extern TLS_CLASS_SPEC integer blas_cpu_number;
void       blas_thread_init(void);

void time_Gemm_nn(
               integer variant, integer type, integer nrepeats, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    m_input, k_input, n_input,
    m, n, k,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    n_repeats,
    variant,
    n_threads,
    n_thread_experiments,
    i, j;

  int n_threads_exp[64];

  char *colors = "brkgmckkk";
  char *ticks =  "o+*xso+*x";
  char m_dim_desc[14];
  char k_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[5];
  char k_dim_tag[5];
  char n_dim_tag[5];
  char nth_str[32];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff,
    d_n;

  FLA_Obj
    A, B, C, C_ref;

  
  /* Initialize FLAME */
  FLA_Init( );


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m k n (-1 means bind to problem size: ", '%' );
  scanf( "%d%d%d", &m_input, &k_input, &n_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, k_input, n_input );
  
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


  if     ( m_input >  0 ) {
    sprintf( m_dim_desc, "m = %d", m_input );
    sprintf( m_dim_tag,  "m%dc", m_input);
  }
  else if( m_input <  -1 ) {
    sprintf( m_dim_desc, "m = p/%d", -m_input );
    sprintf( m_dim_tag,  "m%dp", -m_input );
  }
  else if( m_input == -1 ) {
    sprintf( m_dim_desc, "m = p" );
    sprintf( m_dim_tag,  "m%dp", 1 );
  }
  if     ( k_input >  0 ) {
    sprintf( k_dim_desc, "k = %d", k_input );
    sprintf( k_dim_tag,  "k%dc", k_input);
  }
  else if( k_input <  -1 ) {
    sprintf( k_dim_desc, "k = p/%d", -k_input );
    sprintf( k_dim_tag,  "k%dp", -k_input );
  }
  else if( k_input == -1 ) {
    sprintf( k_dim_desc, "k = p" );
    sprintf( k_dim_tag,  "k%dp", 1 );
  }
  if     ( n_input >  0 ) {
    sprintf( n_dim_desc, "n = %d", n_input );
    sprintf( n_dim_tag,  "n%dc", n_input);
  }
  else if( n_input <  -1 ) {
    sprintf( n_dim_desc, "n = p/%d", -n_input );
    sprintf( n_dim_tag,  "n%dp", -n_input );
  }
  else if( n_input == -1 ) {
    sprintf( n_dim_desc, "n = p" );
    sprintf( n_dim_tag,  "n%dp", 1 );
  }

  m = p_last;
  k = p_last;
  n = p_last;
  



  sprintf( nth_str, "OMP_NUM_THREADS=%d", n_threads_exp[ n_thread_experiments-1 ] );
  putenv( nth_str );
  blas_cpu_number = n_threads_exp[ n_thread_experiments-1 ];
  blas_thread_init();



  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( k < 0 ) k = p / f2c_abs(k_input);
    if( n < 0 ) n = p / f2c_abs(n_input);
	
    FLA_Obj_create( FLA_DOUBLE, m, k, &A );
    FLA_Obj_create( FLA_DOUBLE, k, n, &B );
    FLA_Obj_create( FLA_DOUBLE, m, n, &C );
    FLA_Obj_create( FLA_DOUBLE, m, n, &C_ref );

	
    /* Generate random matrices A, C */
	if( p > 4000 ){
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );
    FLA_Random_matrix( C );
	
    FLA_Copy_external( C, C_ref );
	}
	


    blas_cpu_number = 1;

    //time_Gemm_nn( 0, FLA_ALG_REFERENCE, n_repeats, p, nb_alg,
    //                A, B, C, C_ref, &dtime, &diff, &gflops );

    //fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    //fflush( stdout );

    for ( j = 0; j < n_thread_experiments; j++ ){

      n_threads = n_threads_exp[j];
      blas_cpu_number = n_threads;

      fprintf( stdout, "data_nth%d( %d, 1:3 ) = [ %d  ", n_threads, i, p );
      fflush( stdout );

      time_Gemm_nn( 0, FLA_ALG_REFERENCE, n_repeats, p, nb_alg,
                    A, B, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &C_ref );



  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  fprintf( stdout, "figure;\n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on;\n" );

  /* Plot the data for the other numbers of threads */
  for ( i = 0; i < n_thread_experiments; i++ ){
    fprintf( stdout, "plot( data_nth%d( :,1 ), data_nth%d( :, 2 ), '%c:%c' ); \n", 
             n_threads_exp[ i ], n_threads_exp[ i ], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_thread_experiments-1; i++ )
    fprintf( stdout, "'%d threads', ... \n", n_threads_exp[ i ] );

  fprintf( stdout, "'%d threads', 'Location', 'Best' ); \n", n_threads_exp[ n_thread_experiments-1 ] );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, n_threads_exp[n_thread_experiments-1] * max_gflops );
  fprintf( stdout, "title( 'Goto BLAS dgemm performance (%s, %s, %s)' );\n", 
           m_dim_desc, k_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc gemm_nn_goto_p_%s_%s_%s.eps\n", m_dim_tag, k_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

