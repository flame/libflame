
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

void time_Apply_QUD_UT(
                 int n_repeats, int mB, int mC, int mD, int n, int n_rhs, int b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops );

                         //R_BC, R_BD, C, D, T, bR_BC, bR_BD, bC, bD, &dtime, &diff, &gflops );

int main(int argc, char *argv[])
{
  int 
    datatype,
    n_input,
    n_rhs_input,
    mB_input, mC_input, mD_input,
    n_rhs,
    mB, mC, mD, n,
    p_first, p_last, p_inc,
    p,
    b_alg,
    variant,
    n_repeats,
    i,
    n_variants = 1;
  
  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    R_BD, R_BC, B, C, D, T, W,
    bR_BD, bR_BC, bB, bC, bD;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter n (-1 means bind to problem size): ", '%' );
  scanf( "%d", &n_input );
  fprintf( stdout, "%c %d\n", '%', n_input );

  fprintf( stdout, "%c enter mB mC mD (-1 means bind to problem size): ", '%' );
  scanf( "%d %d %d", &mB_input, &mC_input, &mD_input );
  fprintf( stdout, "%c %d %d %d\n", '%', mB_input, mC_input, mD_input );

  fprintf( stdout, "%c enter n_rhs (-1 means bind to problem size): ", '%' );
  scanf( "%d", &n_rhs_input );
  fprintf( stdout, "%c %d\n", '%', n_rhs_input );


  fprintf( stdout, "\nclear all;\n\n" );



  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    mB    = mB_input;
    mC    = mC_input;
    mD    = mD_input;
    n     = n_input;
    n_rhs = n_rhs_input;

    if( mB    < 0 ) mB    = p / abs(mB_input);
    if( mC    < 0 ) mC    = p / abs(mC_input);
    if( mD    < 0 ) mD    = p / abs(mD_input);
    if( n     < 0 ) n     = p / abs(n_input);
    if( n_rhs < 0 ) n_rhs = p / abs(n_rhs_input);

    for ( variant = 0; variant < n_variants; variant++ ){
      
      FLA_Obj_create( datatype, mB, n, 0, 0, &B );
      FLA_Obj_create( datatype, mC, n, 0, 0, &C );
      FLA_Obj_create( datatype, mD, n, 0, 0, &D );
      FLA_Obj_create( datatype, b_alg, n, 0, 0, &T );
      FLA_Obj_create( datatype, n,  n, 0, 0, &R_BC );
      FLA_Obj_create( datatype, n,  n, 0, 0, &R_BD );

      FLA_Obj_create( datatype, mB, n_rhs, 0, 0, &bB );
      FLA_Obj_create( datatype, mC, n_rhs, 0, 0, &bC );
      FLA_Obj_create( datatype, mD, n_rhs, 0, 0, &bD );
      FLA_Obj_create( datatype, n,  n_rhs, 0, 0, &bR_BC );
      FLA_Obj_create( datatype, n,  n_rhs, 0, 0, &bR_BD );

      FLA_Apply_QUD_UT_create_workspace( T, bR_BD, &W );

      FLA_Random_matrix( B );
      FLA_Random_matrix( C );
      FLA_Random_matrix( D );

      FLA_Random_matrix( bB );
      FLA_Random_matrix( bC );
      FLA_Random_matrix( bD );

      FLA_Set( FLA_ZERO, R_BD );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R_BD );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D, FLA_ONE, R_BD );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R_BD );

      FLA_Set( FLA_ZERO, R_BC );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R_BC );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C, FLA_ONE, R_BC );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R_BC );

      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BD );
      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, D, bD, FLA_ONE,  bR_BD );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BD, bR_BD );

      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BC );
      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, C, bC, FLA_ONE,  bR_BC );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BC, bR_BC );

      FLA_UDdate_UT( R_BD, C, D, T );

      fprintf( stdout, "data_apqud_ut( %d, 1:5 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Apply_QUD_UT( n_repeats, mB, mC, mD, n, n_rhs, b_alg,
                         R_BC, R_BD, C, D, T, W, bR_BC, bR_BD, bC, bD, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &B );
      FLA_Obj_free( &C );
      FLA_Obj_free( &D );
      FLA_Obj_free( &T );
      FLA_Obj_free( &W );
      FLA_Obj_free( &R_BC );
      FLA_Obj_free( &R_BD );

      FLA_Obj_free( &bB );
      FLA_Obj_free( &bC );
      FLA_Obj_free( &bD );
      FLA_Obj_free( &bR_BC );
      FLA_Obj_free( &bR_BD );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_variants; i++ ) {
    fprintf( stdout, "plot( data_apqudut( :,1 ), data_apqudut( :, 2 ), '%c:%c' ); \n",
            colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_apqudut( :,1 ), data_apqudut( :, 4 ), '%c-.%c' ); \n",
            colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_variants; i++ )
    fprintf( stdout, "'ref\\_qr\\_ut', 'fla\\_qr\\_ut', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Apply_QUD_UT front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc apqudut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

