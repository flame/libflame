/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

void time_UDdate_UT_inc(
                 int variant, int type, int n_repeats, int mB, int mC, int mD, int n,
                 FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_input,
    mB_input, mC_input, mD_input,
    mB, mC, mD, n,
    p_first, p_last, p_inc,
    p,
    variant,
    n_repeats,
    i,
    n_variants = 1;

  dim_t
    b_flash,
    b_alg,
    n_threads;
  
  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    B_flat, C_flat, D_flat, R_flat, E_flat,
    C, D, T, W, R, E;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize:", '%' );
  scanf( "%u", &b_alg );
  fprintf( stdout, "%c %u\n", '%', b_alg );

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
    mB = mB_input;
    mC = mC_input;
    mD = mD_input;
    n  = n_input;

    if( mB < 0 ) mB = p / abs(mB_input);
    if( mC < 0 ) mC = p / abs(mC_input);
    if( mD < 0 ) mD = p / abs(mD_input);
    if( n  < 0 ) n  = p / abs(n_input);

    for ( variant = 0; variant < n_variants; variant++ ){
      
      FLA_Obj_create( datatype, mB, n, 0, 0, &B_flat );
      FLA_Obj_create( datatype, mC, n, 0, 0, &C_flat );
      FLA_Obj_create( datatype, mD, n, 0, 0, &D_flat );

      FLA_Random_matrix( B_flat );
      FLA_Random_matrix( C_flat );
      FLA_Random_matrix( D_flat );

      FLA_Obj_create( datatype, n, n, 0, 0, &R_flat );
      FLA_Obj_create( datatype, n, n, 0, 0, &E_flat );

      FLA_Set( FLA_ZERO, R_flat );
      FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B_flat, FLA_ONE, R_flat );
      FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D_flat, FLA_ONE, R_flat );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R_flat );

      FLA_Set( FLA_ZERO, E_flat );
      FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B_flat, FLA_ONE, E_flat );
      FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C_flat, FLA_ONE, E_flat );
      FLA_Chol( FLA_UPPER_TRIANGULAR, E_flat );

      FLASH_UDdate_UT_inc_create_hier_matrices( R_flat, C_flat, D_flat, 1, &b_flash, b_alg,
                                                &R, &C, &D, &T, &W );
      FLASH_Obj_create_hier_copy_of_flat( E_flat, 1, &b_flash, &E );


      fprintf( stdout, "data_uddate_ut_inc( %d, 1:5 ) = [ %d  ", i, p );
      fflush( stdout );

      time_UDdate_UT_inc( variant, FLA_ALG_FRONT, n_repeats, mB, mC, mD, n,
                          C, D, T, W, R, E, &dtime, &diff, &gflops );


      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &B_flat );
      FLASH_Obj_free( &C_flat );
      FLASH_Obj_free( &D_flat );
      FLASH_Obj_free( &R_flat );
      FLASH_Obj_free( &E_flat );
      FLASH_Obj_free( &C );
      FLASH_Obj_free( &D );
      FLASH_Obj_free( &T );
      FLASH_Obj_free( &W );
      FLASH_Obj_free( &R );
      FLASH_Obj_free( &E );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_variants; i++ ) {
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 2 ), '%c:%c' ); \n",
            colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 4 ), '%c-.%c' ); \n",
            colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_variants; i++ )
    fprintf( stdout, "'ref\\_qr\\_ut', 'fla\\_qr\\_ut', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME UDdate_UT front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc qr_ut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

