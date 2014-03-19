/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_GEMM_FULL     0
#define FLA_ALG_GEMM_FLOPEQ   1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_UNB_ASM       3
#define FLA_ALG_BLOCKED       4

void time_Apply_G_rf(
               int variant, int type, int n_repeats, int m, int k, int n, int b_alg,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj G,
               FLA_Obj Ag, FLA_Obj Bg, FLA_Obj Cg, FLA_Obj Dg, FLA_Obj Eg,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, k_input, n_input,
    m, k, n,
    p_first, p_last, p_inc,
    p,
    b_alg,
    n_repeats,
    i, mr, nr, kr, kf,
    datatype, dt_real, dt_comp;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, A_ref, G, Ag, Bg, Cg, Dg, Eg;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n k (-1 means bind to problem size): ", '%' );
  scanf( "%d %d %d", &m_input, &n_input, &k_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, n_input, k_input );


  fprintf( stdout, "\n" );



  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( k < 0 ) k = p / abs(k_input);
    if( n < 0 ) n = p / abs(n_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;


    // For apply_g
    FLA_Obj_create( datatype, m,   n, 0, 0, &A );
    FLA_Obj_create( datatype, m,   n, 0, 0, &A_ref );

	if ( FLA_Obj_is_double_precision( A ) ) dt_comp = FLA_DOUBLE_COMPLEX;
	else                                    dt_comp = FLA_COMPLEX;

    FLA_Obj_create( dt_comp, n-1, k, 0, 0, &G );

    // For gemm
    dt_real = FLA_Obj_datatype_proj_to_real( A );
    mr      =     m;
    nr      =     n;
    kr      = 3 * k;
    kf      =     n;

    FLA_Obj_create( dt_real,  mr, nr, 0, 0, &Ag );

    FLA_Obj_create( dt_real,  mr, kf, 0, 0, &Bg );
    FLA_Obj_create( dt_real,  kf, nr, 0, 0, &Cg );

    FLA_Obj_create( dt_real,  mr, kr, 0, 0, &Dg );
    FLA_Obj_create( dt_real,  kr, nr, 0, 0, &Eg );

    FLA_Random_matrix( A );
    //FLA_Set_to_identity( A );
    //FLA_Set( FLA_ZERO, A );
    //FLA_Set_diag( FLA_TWO, A );
    //FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    //FLA_Random_tri_matrix( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    FLA_Random_matrix( G );
    //fill_cs( G );

    FLA_Random_matrix( Ag );
    FLA_Random_matrix( Bg );
    FLA_Random_matrix( Cg );
    FLA_Random_matrix( Dg );
    FLA_Random_matrix( Eg );

    time_Apply_G_rf( 0, FLA_ALG_GEMM_FULL, n_repeats, mr, kf, nr, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_gemm_full(     %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 0, FLA_ALG_GEMM_FLOPEQ, n_repeats, mr, kr, nr, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_gemm_rank3k(   %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 1, FLA_ALG_UNB_OPT, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_opt_var1( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 1, FLA_ALG_UNB_ASM, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_asm_var1( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 2, FLA_ALG_BLOCKED, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_blk_var2( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 6, FLA_ALG_BLOCKED, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_blk_var6( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 9, FLA_ALG_BLOCKED, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_blk_var9( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Apply_G_rf( 3, FLA_ALG_BLOCKED, n_repeats, m, k, n, b_alg,
                     A, A_ref, G, Ag, Bg, Cg, Dg, Eg, &dtime, &diff, &gflops );

    fprintf( stdout, "data_appg_blk_var3( %d, 1:3 ) = [ %4d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );



    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &A_ref );
    FLA_Obj_free( &G );
    FLA_Obj_free( &Ag );
    FLA_Obj_free( &Bg );
    FLA_Obj_free( &Cg );
    FLA_Obj_free( &Dg );
    FLA_Obj_free( &Eg );
  }

  FLA_Finalize( );

  return 0;
}

