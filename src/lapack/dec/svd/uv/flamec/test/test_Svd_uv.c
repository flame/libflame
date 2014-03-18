
#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3


void time_Svd_uv(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj U, FLA_Obj V, FLA_Obj s,
               double *dtime, double *diff1, double* diff2, double *gflops, int* k_perf );


int main(int argc, char *argv[])
{
  int 
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    k_accum,
    b_alg,
    n_iter_max,
    min_m_n,
    variant,
    n_repeats,
    i,
    k_perf,
    n_variants = 2;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];
  char n_dim_desc[14];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff1,
    diff2;

  FLA_Datatype datatype, dt_real;

  FLA_Obj
    A, s, U, V, alpha;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter n_iter_max:", '%' );
  scanf( "%d", &n_iter_max );
  fprintf( stdout, "%c %d\n", '%', n_iter_max );

  fprintf( stdout, "%c enter number of Givens rotations to accumulate:", '%' );
  scanf( "%d", &k_accum );
  fprintf( stdout, "%c %d\n", '%', k_accum );

  fprintf( stdout, "%c enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


  fprintf( stdout, "\n" );


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
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( n < 0 ) n = p / abs(n_input);

    min_m_n = min( m, n );

    //datatype = FLA_FLOAT;
    datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    //datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,       n, 0, 0, &A );
    FLA_Obj_create( datatype, m,       m, 0, 0, &U );
    FLA_Obj_create( datatype, n,       n, 0, 0, &V );

    dt_real = FLA_Obj_datatype_proj_to_real( A );

    FLA_Obj_create( dt_real,  min_m_n, 1, 0, 0, &s );
    FLA_Obj_create( dt_real,  1,       1, 0, 0, &alpha );

    FLA_Random_unitary_matrix( U );
    FLA_Random_unitary_matrix( V );

    FLA_Fill_with_linear_dist( FLA_ZERO, FLA_ONE, s );

    //FLA_Fill_with_inverse_dist( FLA_ONE, s );

    //*FLA_DOUBLE_PTR( alpha ) = 1.0 / sqrt( (double) min_m_n );
    //FLA_Fill_with_geometric_dist( alpha,   s );

    {
      FLA_Obj UL, UR;
      FLA_Obj VL, VR;

      FLA_Part_1x2( U,   &UL, &UR,   min_m_n, FLA_LEFT );
      FLA_Part_1x2( V,   &VL, &VR,   min_m_n, FLA_LEFT );

      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, UL, VL, FLA_ZERO, A );
    }

/*
*FLA_DOUBLE_PTR( alpha ) = 1.0e-169;
FLA_Scal( alpha, A );
FLA_Obj_show( "A", A, "%10.2e", "" );
*/

    FLA_Set( FLA_ZERO, s );
    FLA_Set( FLA_ZERO, U );
    FLA_Set( FLA_ZERO, V );

    time_Svd_uv( 0, FLA_ALG_REFERENCE, n_repeats, m, n, n_iter_max, k_accum, b_alg,
                  A, s, U, V, &dtime, &diff1, &diff2, &gflops, &k_perf );

    fprintf( stdout, "data_REFq( %d, 1:5 ) = [ %d %6.3lf %8.2e  %6.2le %6.2le %5d ]; \n", i, p, gflops, dtime, diff1, diff2, k_perf );
    fflush( stdout );

    time_Svd_uv( -1, FLA_ALG_REFERENCE, n_repeats, m, n, n_iter_max, k_accum, b_alg,
                  A, s, U, V, &dtime, &diff1, &diff2, &gflops, &k_perf );

    fprintf( stdout, "data_REFd( %d, 1:5 ) = [ %d %6.3lf %8.2e  %6.2le %6.2le %5d ]; \n", i, p, gflops, dtime, diff1, diff2, k_perf );
    fflush( stdout );

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:5 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Svd_uv( variant, FLA_ALG_UNBLOCKED, n_repeats, m, n, n_iter_max, k_accum, b_alg,
                    A, s, U, V, &dtime, &diff1, &diff2, &gflops, &k_perf );

      fprintf( stdout, "%6.3lf %8.2e  %6.2le %6.2le %5d ", gflops, dtime, diff1, diff2, k_perf );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &s );
    FLA_Obj_free( &U );
    FLA_Obj_free( &V );
    FLA_Obj_free( &alpha );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 1; i <= n_variants; i++ ) {
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i < n_variants; i++ )
    fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', ... \n", i, i );
  fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d' ); \n", i, i );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Svd_uv performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

