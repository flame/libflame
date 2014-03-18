#include <time.h>

#include "FLAME.h"

#define N_VARIANTS 5

#define FLA_ALG_REFERENCE 0


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
    n_repeats,
    variant,
    i, j;

  int n_threads_exp[64];

  char *colors = "brkgmckkk";
  char *ticks  = "o+*xso+*x";
  char m_dim_desc[14] = "";
  char k_dim_desc[14] = "";
  char n_dim_desc[14] = "";
  char m_dim_tag[5] = "";
  char k_dim_tag[5] = "";
  char n_dim_tag[5] = "";

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff,
    d_n;

  FLA_Obj
    A, B, C, C_ref,
    ATL, ATR, ABL, ABR,
    BTL, BTR, BBL, BBR,
    CTL, CTR, CBL, CBR;
  
  /* Initialize FLAME */
  FLA_Init( );


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m k n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d%d", &m_input, &k_input, &n_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, k_input, n_input );
  
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

  FLA_Obj_create( FLA_DOUBLE, m, k, &A );
  FLA_Obj_create( FLA_DOUBLE, k, n, &B );
  FLA_Obj_create( FLA_DOUBLE, m, n, &C );
  FLA_Obj_create( FLA_DOUBLE, m, n, &C_ref );




  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    
    m = m_input;
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( k < 0 ) k = p / abs(k_input);
    if( n < 0 ) n = p / abs(n_input);


    FLA_Part_2x2( A, &ATL, /**/ &ATR,
                  /* *************** */
                     &ABL, /**/ &ABR,
               m, k, FLA_TL );

    FLA_Part_2x2( B, &BTL, /**/ &BTR,
                  /* *************** */
                     &BBL, /**/ &BBR,
               k, n, FLA_TL );

    FLA_Part_2x2( C, &CTL, /**/ &CTR,
                  /* *************** */
                     &CBL, /**/ &CBR,
               m, n, FLA_TL );


    FLA_Random_matrix( ATL );
    FLA_Random_matrix( BTL );
    FLA_Random_matrix( CTL );


    time_Gemm_nn( 0, FLA_ALG_REFERENCE, n_repeats, n, nb_alg,
		ATL, BTL, CTL, C_ref, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );

    fprintf( stdout, "\n" );


  }


  FLA_Obj_free( &A );
  FLA_Obj_free( &B );
  FLA_Obj_free( &C );
  FLA_Obj_free( &C_ref );

  FLA_Finalize( );

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), 'b-' ); \n" );

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Goto BLAS gemm', 2 ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n");
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'Goto BLAS gemm\\_nn performance (%s, %s, %s)' );\n", 
           m_dim_desc, k_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc gemm_nn_goto_%s_%s_%s.eps\n", m_dim_tag, k_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
}

