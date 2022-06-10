/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

void time_Apply_QUD_UT_inc(
                 int n_repeats, int mB, int mC, int mD, int n, int n_rhs, dim_t b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops );


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
    variant,
    n_repeats,
    i,
    n_variants = 1;
  
  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  dim_t b_alg, b_flash, n_threads;

  FLA_Obj
    R_BD_flat, R_BC_flat, B_flat, C_flat, D_flat,
    R_BD, R_BC, B, C, D, T, W, W2,
    bR_BD, bR_BC, bB, bC, bD;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize:", '%' );
  scanf( "%lu", &b_alg );
  fprintf( stdout, "%c %lu\n", '%', b_alg );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &b_flash );
  fprintf( stdout, "%c %u\n", '%', b_flash );

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

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%u", &n_threads );
  fprintf( stdout, "%c %u\n", '%', n_threads );


  fprintf( stdout, "\nclear all;\n\n" );



  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );
  //FLASH_Queue_set_verbose_output( TRUE );
  //FLASH_Queue_disable();

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    mB    = mB_input;
    mC    = mC_input;
    mD    = mD_input;
    n     = n_input;
    n_rhs = n_rhs_input;

    if( mB    < 0 ) mB    = p / f2c_abs(mB_input);
    if( mC    < 0 ) mC    = p / f2c_abs(mC_input);
    if( mD    < 0 ) mD    = p / f2c_abs(mD_input);
    if( n     < 0 ) n     = p / f2c_abs(n_input);
    if( n_rhs < 0 ) n_rhs = p / f2c_abs(n_rhs_input);

    for ( variant = 0; variant < n_variants; variant++ ){
      
      FLA_Obj_create( datatype, mB, n, 0, 0, &B_flat );
      FLA_Obj_create( datatype, mC, n, 0, 0, &C_flat );
      FLA_Obj_create( datatype, mD, n, 0, 0, &D_flat );
      FLA_Obj_create( datatype, n,  n, 0, 0, &R_BC_flat );
      FLA_Obj_create( datatype, n,  n, 0, 0, &R_BD_flat );

      FLA_Random_matrix( B_flat );
      FLA_Random_matrix( C_flat );
      FLA_Random_matrix( D_flat );

      FLA_Set( FLA_ZERO, R_BD_flat );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_ONE, B_flat, FLA_ONE, R_BD_flat );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_ONE, D_flat, FLA_ONE, R_BD_flat );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R_BD_flat );

      FLA_Set( FLA_ZERO, R_BC_flat );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_ONE, B_flat, FLA_ONE, R_BC_flat );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                         FLA_ONE, C_flat, FLA_ONE, R_BC_flat );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R_BC_flat );

      FLASH_Obj_create_hier_copy_of_flat( B_flat, 1, &b_flash, &B );
      FLASH_Obj_create_hier_copy_of_flat( R_BC_flat, 1, &b_flash, &R_BC );
      FLASH_UDdate_UT_inc_create_hier_matrices( R_BD_flat, C_flat, D_flat,
                                                1, &b_flash, b_alg, &R_BD, &C, &D, &T, &W );

      FLASH_Obj_create( datatype, mB, n_rhs, 1, &b_flash, &bB );
      FLASH_Obj_create( datatype, mC, n_rhs, 1, &b_flash, &bC );
      FLASH_Obj_create( datatype, mD, n_rhs, 1, &b_flash, &bD );
      FLASH_Obj_create( datatype, n,  n_rhs, 1, &b_flash, &bR_BC );
      FLASH_Obj_create( datatype, n,  n_rhs, 1, &b_flash, &bR_BD );

      FLASH_Random_matrix( bB );
      FLASH_Random_matrix( bC );
      FLASH_Random_matrix( bD );

      FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BD );
      FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, D, bD, FLA_ONE,  bR_BD );
      FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BD, bR_BD );

      FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BC );
      FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, C, bC, FLA_ONE,  bR_BC );
      FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BC, bR_BC );

      FLASH_UDdate_UT_inc( R_BD, C, D, T, W );

      FLASH_Apply_QUD_UT_inc_create_workspace( T, bR_BD, &W2 );

      fprintf( stdout, "data_apqudutinc( %d, 1:5 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Apply_QUD_UT_inc( n_repeats, mB, mC, mD, n, n_rhs, b_alg,
                             R_BC, R_BD, C, D, T, W2, bR_BC, bR_BD, bC, bD, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &B_flat );
      FLA_Obj_free( &C_flat );
      FLA_Obj_free( &D_flat );
      FLA_Obj_free( &R_BC_flat );
      FLA_Obj_free( &R_BD_flat );

      FLASH_Obj_free( &B );
      FLASH_Obj_free( &C );
      FLASH_Obj_free( &D );
      FLASH_Obj_free( &T );
      FLASH_Obj_free( &W );
      FLASH_Obj_free( &W2 );
      FLASH_Obj_free( &R_BC );
      FLASH_Obj_free( &R_BD );

      FLASH_Obj_free( &bB );
      FLASH_Obj_free( &bC );
      FLASH_Obj_free( &bD );
      FLASH_Obj_free( &bR_BC );
      FLASH_Obj_free( &bR_BD );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_variants; i++ ) {
    fprintf( stdout, "plot( data_apqudutinc( :,1 ), data_apqudutinc( :, 2 ), '%c:%c' ); \n",
            colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_apqudutinc( :,1 ), data_apqudutinc( :, 4 ), '%c-.%c' ); \n",
            colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_variants; i++ )
    fprintf( stdout, "'ref\\_qr\\_ut', 'fla\\_qr\\_ut', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Apply_QUD_UT_inc front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc apqudut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

