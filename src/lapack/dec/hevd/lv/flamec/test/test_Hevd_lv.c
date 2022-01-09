/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3


void time_Hevd_lv(
               integer variant, integer type, integer n_repeats, integer m, integer n_iter_max, integer k_accum, integer b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops, integer* k_perf );


int main(int argc, char *argv[])
{
  integer 
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    k_accum,
    b_alg,
    n_iter_max,
    variant,
    n_repeats,
    i,
    k_perf,
    n_variants = 2;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff1,
    diff2;

  FLA_Datatype datatype, dt_real;

  FLA_Obj
    A, l, l_save, Q, Ql, alpha, n_clust;
  

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

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );


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


  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    //datatype = FLA_FLOAT;
    datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    //datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,  m, 0, 0, &A );
    FLA_Obj_create( datatype, m,  m, 0, 0, &Q );
    FLA_Obj_create( datatype, m,  m, 0, 0, &Ql );

    dt_real = FLA_Obj_datatype_proj_to_real( A );

    FLA_Obj_create( dt_real,  m,  1, 0, 0, &l );
    FLA_Obj_create( dt_real,  m,  1, 0, 0, &l_save );
    FLA_Obj_create( dt_real,  1,  1, 0, 0, &alpha );
    FLA_Obj_create( FLA_INT,  1,  1, 0, 0, &n_clust );

    FLA_Random_unitary_matrix( Q );

	FLA_Fill_with_linear_dist( FLA_ZERO, FLA_ONE, l );

	//FLA_Fill_with_inverse_dist( FLA_ONE, l );

    //*FLA_DOUBLE_PTR( alpha ) = 1.0 / sqrt( (double) m );
    //*FLA_DOUBLE_PTR( alpha ) = 0.1;
	//FLA_Fill_with_geometric_dist( alpha,   l );

    //*FLA_DOUBLE_PTR( alpha ) = 1.2;
	//FLA_Fill_with_logarithmic_dist( alpha,   l );

    //*FLA_DOUBLE_PTR( alpha ) = (double)m;
	//FLA_Fill_with_random_dist( FLA_ZERO, alpha, l );

    //*FLA_INT_PTR( n_clust ) =  8;
    //*FLA_DOUBLE_PTR( alpha ) = (double) 1.0e-03;
	//FLA_Fill_with_cluster_dist( n_clust, alpha, l );

    FLA_Copy( l, l_save );

    {
      FLA_Copy( Q, Ql );
      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, Ql );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, Ql, Q, FLA_ZERO, A );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    }


    FLA_Set( FLA_ZERO, l );


    time_Hevd_lv( 0, FLA_ALG_REFERENCE, n_repeats, m, n_iter_max, k_accum, b_alg,
                  A, l, &dtime, &diff1, &diff2, &gflops, &k_perf );

    fprintf( stdout, "data_REFq( %d, 1:3 ) = [ %d %6.3lf %8.2e  %6.2le %6.2le %5d ]; \n", i, p, gflops, dtime, diff1, diff2, k_perf );
    fflush( stdout );

    time_Hevd_lv( -1, FLA_ALG_REFERENCE, n_repeats, m, n_iter_max, k_accum, b_alg,
                  A, l, &dtime, &diff1, &diff2, &gflops, &k_perf );

    fprintf( stdout, "data_REFd( %d, 1:3 ) = [ %d %6.3lf %8.2e  %6.2le %6.2le %5d ]; \n", i, p, gflops, dtime, diff1, diff2, k_perf );
    fflush( stdout );

    time_Hevd_lv( -2, FLA_ALG_REFERENCE, n_repeats, m, n_iter_max, k_accum, b_alg,
                  A, l, &dtime, &diff1, &diff2, &gflops, &k_perf );

    fprintf( stdout, "data_REFr( %d, 1:3 ) = [ %d %6.3lf %8.2e  %6.2le %6.2le %5d ]; \n", i, p, gflops, dtime, diff1, diff2, k_perf );
    fflush( stdout );


    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:3 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Hevd_lv( variant, FLA_ALG_UNBLOCKED, n_repeats, m, n_iter_max, k_accum, b_alg,
                    A, l, &dtime, &diff1, &diff2, &gflops, &k_perf );

/*
if ( m == 145 )
{
printf( "max elemwise diff = %e\n", FLA_Max_elemwise_diff( l_save, l ) );
FLA_Axpy( FLA_MINUS_ONE, l_save, l );
FLA_Obj_show( "eigval diff:", l, "%17.10e", "" );
}
*/

      fprintf( stdout, "%6.3lf %8.2e  %6.2le %6.2le %5d ", gflops, dtime, diff1, diff2, k_perf );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &Q );
    FLA_Obj_free( &Ql );
    FLA_Obj_free( &l );
    FLA_Obj_free( &l_save );
    FLA_Obj_free( &alpha );
    FLA_Obj_free( &n_clust );
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
  fprintf( stdout, "title( 'FLAME Hevd_lv performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

