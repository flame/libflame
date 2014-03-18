#include <time.h>

#include "FLAME.h"

#define FLA_ALG_REFERENCE    0
#define FLA_ALG_OPENMP_1TASK     5
#define FLA_ALG_OPENMP_2TASKS    6
#define FLA_ALG_OPENMP_2LOOPS    7



void time_Syrk_ln(
               int variant, int type, int n_repeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, k_input, n_input,
    m, n, k,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    variant,
    n_repeats,
    n_thread_experiments,
    i, j;
  
  int  n_threads[64];
  int  blocksize_tag0[64];
  int  blocksize_tag1[64];

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char k_dim_desc[14];
  char m_dim_tag[9];
  char k_dim_tag[9];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, B, C, C_ref;
  
  /* Initialize FLAME */
  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m k (-1 means bind to problem size): ", '%' );
  scanf( "%d%d", &m_input, &k_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, k_input );

  fprintf( stdout, "%c enter number of thread experiments: ", '%' );
  scanf( "%d", &n_thread_experiments );
  fprintf( stdout, "%c %d\n", '%', n_thread_experiments );

  fprintf( stdout, "%c enter t for each experiment: ", '%' );
  for( i = 0; i < n_thread_experiments; ++i )
    scanf( "%d", &n_threads[i] );

  fprintf( stdout, "\n" );
  for( i = 0; i < n_thread_experiments; ++i )
    fprintf( stdout, "%c %2d threads\n", '%', n_threads[i] );



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




  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    k = k_input;

    if( m < 0 ) m = p / abs(m_input);
    if( k < 0 ) k = p / abs(k_input);


    /* Allocate space for the matrices */
    FLA_Obj_create( FLA_DOUBLE, m, k, &A );
    FLA_Obj_create( FLA_DOUBLE, m, m, &C );
    FLA_Obj_create( FLA_DOUBLE, m, m, &C_ref );

    /* Generate random matrices A, C */
    FLA_Random_matrix( A );
    FLA_Random_matrix( C );

    FLA_Copy_external( C, C_ref );


    /* Time the reference implementation */
    time_Syrk_ln( 0, FLA_ALG_REFERENCE, n_repeats, n, nb_alg,
                  A, B, C, C_ref, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );

    FLA_Obj_free( &A );
    FLA_Obj_free( &C );
    FLA_Obj_free( &C_ref );
    fprintf( stdout, "\n" );


  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  fprintf( stdout, "figure;\n" );

  /* Plot the performance of the reference implementation */
  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on;\n" );

  fprintf( stdout, "legend( Reference ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, n_threads[n_thread_experiments-1] * max_gflops );
  fprintf( stdout, "title( 'GotoBLAS syrk\\_ln performance (%s, %s)' );\n", 
           m_dim_desc, k_dim_desc );
  fprintf( stdout, "print -depsc syrk_ln_goto_p_%s_%s.eps\n", m_dim_tag, k_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

