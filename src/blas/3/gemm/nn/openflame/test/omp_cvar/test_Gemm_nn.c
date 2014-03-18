#include <time.h>

#include "FLAME.h"
#include "FLA_task_partitioning.h"

#define N_VARIANTS 5

#define FLA_ALG_REFERENCE    0
#define FLA_ALG_OPENMP_BVAR  1
#define FLA_ALG_OPENMP_CVAR  2


void time_Gemm_nn(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, k_input, n_input,
    m, n, k,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    nrepeats,
    variant,
    n_threads,
    n_thread_experiments,
    i, j,
    nvariants = N_VARIANTS;
  
  int  n_threads_exp[64];
  int  n_threads_exp_m[64];
  int  n_threads_exp_k[64];
  int  n_threads_exp_n[64];

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char k_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[5];
  char k_dim_tag[5];
  char n_dim_tag[5];

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
  FLA_Task_partitioning_init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &nrepeats );
  fprintf( stdout, "%c %d\n", '%', nrepeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m k n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d%d", &m_input, &k_input, &n_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, k_input, n_input );

  fprintf( stdout, "%c enter variant or variant-permutation (1..6,13,31,15,35): ", '%' );
  scanf( "%d", &variant );
  fprintf( stdout, "%c %d\n", '%', variant );

  fprintf( stdout, "%c enter number of thread experiments: ", '%' );
  scanf( "%d", &n_thread_experiments );
  fprintf( stdout, "%c %d\n", '%', n_thread_experiments );

  fprintf( stdout, "%c enter t, t_m, t_k, and t_n for each experiment: ", '%' );
  for( i = 0; i < n_thread_experiments; ++i )
    scanf( "%d %d %d %d", &n_threads_exp[i], &n_threads_exp_m[i], &n_threads_exp_k[i], &n_threads_exp_n[i] );

  fprintf( stdout, "\n" );
  for( i = 0; i < n_thread_experiments; ++i )
    fprintf( stdout, "%c %2d = %2d x %2d x %2d\n", '%', n_threads_exp[i], n_threads_exp_m[i], n_threads_exp_k[i], n_threads_exp_n[i] );



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




  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( k < 0 ) k = p / abs(k_input);
    if( n < 0 ) n = p / abs(n_input);


    /* Allocate space for the matrices */
    FLA_Obj_create( FLA_DOUBLE, m, k, &A );
    FLA_Obj_create( FLA_DOUBLE, k, n, &B );
    FLA_Obj_create( FLA_DOUBLE, m, n, &C );
    FLA_Obj_create( FLA_DOUBLE, m, n, &Cref );

    /* Generate random matrices A, C */
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );
    FLA_Random_matrix( C );

    FLA_Copy_external( C, Cref );


    /* Time the reference implementation */
    time_Gemm_nn( 0, FLA_ALG_REFERENCE, nrepeats, n, nb_alg,
                  A, B, C, Cref, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );

    for ( j = 0; j < n_thread_experiments; j++ ){
      
      n_threads = n_threads_exp[j];
      FLA_Task_partitioning_set( n_threads_exp[j], n_threads_exp_m[j], n_threads_exp_k[j], n_threads_exp_n[j] );
      FLA_omp_set_num_threads( n_threads_exp[j] );
      FLA_omp_set_num_stages( n_threads_exp_k[j] );

      fprintf( stdout, "data_nth%d_%dx%dx%d( %d, 1:3 ) = [ %d  ", 
               n_threads, n_threads_exp_m[j], n_threads_exp_k[j], n_threads_exp_n[j], i, p );
      fflush( stdout );

      //time_Gemm_nn( variant, FLA_ALG_OPENMP_BVAR, nrepeats, n, nb_alg,
      time_Gemm_nn( variant, FLA_ALG_OPENMP_CVAR, nrepeats, p, nb_alg,
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
  //fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on;\n" );

  /* Plot the data for the other numbers of threads */
  for ( i = 0; i < n_thread_experiments; i++ ){
    fprintf( stdout, "plot( data_nth%d_%dx%dx%d( :,1 ), data_nth%d_%dx%dx%d( :, 2 ), '%c:%c' ); \n", 
	    n_threads_exp[ i ], n_threads_exp_m[i], n_threads_exp_k[i], n_threads_exp_n[i],
        n_threads_exp[ i ], n_threads_exp_m[i], n_threads_exp_k[i], n_threads_exp_n[i],
        colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_thread_experiments-1; i++ )
    fprintf( stdout, "'n\\_threads %d=%dx%dx%d', ... \n", 
             n_threads_exp[ i ], n_threads_exp_m[i], n_threads_exp_k[i], n_threads_exp_n[i] );

  i = n_thread_experiments-1;
  fprintf( stdout, "'n\\_threads %d=%dx%dx%d', 2 ); \n", 
             n_threads_exp[ i ], n_threads_exp_m[i], n_threads_exp_k[i], n_threads_exp_n[i] );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, n_threads_exp[n_thread_experiments-1] * max_gflops );
  fprintf( stdout, "title( 'OpenFLAME gemm\\_nn\\_var%d performance (%s, %s, %s)' );\n", 
           variant, m_dim_desc, k_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc gemm_nn_ompfac_var%d_%s_%s_%s.eps\n", variant, m_dim_tag, k_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

