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

FLA_Error REF_Hevd_lv( FLA_Obj A, FLA_Obj l );
FLA_Error REF_Hevdd_lv( FLA_Obj A, FLA_Obj l );
FLA_Error REF_Hevdr_lv( FLA_Obj A, FLA_Obj l, FLA_Obj Z );
void time_Hevd_lv(
               int variant, int type, int n_repeats, int m, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops, int* k_perf );


void time_Hevd_lv(
               int variant, int type, int n_repeats, int m, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops, int* k_perf )
{
  int irep;

  double
    k, dtime_old = 1.0e9;

  FLA_Obj
    A_save, Z;

  if (
       //( variant == 0 ) ||
       //( variant == -1 ) ||
       //( variant == -2 ) ||
       //( variant == 1 ) ||
       //( variant == 2 ) ||
       FALSE
     )
  {
    *gflops = 0.0;
    *dtime  = 0.0;
    *diff1  = 0.0;
    *diff2  = 0.0;
    *k_perf  = 0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  FLA_Copy_external( A, A_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Hevd_lv( A, l );
      *k_perf = 0;
      break;

    case -1:
      REF_Hevdd_lv( A, l );
      *k_perf = 0;
      break;

    case -2:
      REF_Hevdr_lv( A, l, Z );
      *k_perf = 0;
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        *k_perf = FLA_Hevd_lv_unb_var1( n_iter_max, A, l, k_accum, b_alg );
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        *k_perf = FLA_Hevd_lv_unb_var2( n_iter_max, A, l, k_accum, b_alg );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }
  {
    FLA_Obj V, A_rev_evd, norm, eye;
    FLA_Obj ATL, ATR, ABL, ABR;
    FLA_Obj lT, lB;

//FLA_Obj_show( "A_save", A_save, "%9.2e + %9.2e ", "" );
//FLA_Obj_show( "A_evd", A, "%9.2e + %9.2e ", "" );

//if ( variant == 2 && m == 145 )
//FLA_Obj_show( "l_evd", l, "%17.10e ", "" );

    if ( variant == -2 ) FLA_Copy( Z, A );

	FLA_Sort_evd( FLA_BACKWARD, l, A );

/*
if ( variant > 0 && m == 145 ){
FLA_Part_2x2( A,   &ATL, &ATR,
                   &ABL, &ABR,   3, 3, FLA_BR );
FLA_Part_2x1( l,   &lT,
                   &lB,    3, FLA_BOTTOM );

FLA_Obj_show( "V", ABR, "%17.10e + %17.10e ", "" );
FLA_Obj_show( "l", lB, "%17.10e ", "" );
}
*/

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &V ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_rev_evd ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &eye ); 
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, A );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, A, V, FLA_ZERO, A_rev_evd );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A_rev_evd );
 
    FLA_Axpy( FLA_MINUS_ONE, A_save, A_rev_evd );
    FLA_Norm_frob( A_rev_evd, norm );
    FLA_Obj_extract_real_scalar( norm, diff1 );

    FLA_Set_to_identity( eye );
	FLA_Copy( V, A_rev_evd );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, V, A_rev_evd, FLA_MINUS_ONE, eye );
    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff2 );

    FLA_Obj_free( &V );
    FLA_Obj_free( &A_rev_evd );
    FLA_Obj_free( &eye );
    FLA_Obj_free( &norm );
  }

  k = 2.00;

  if ( FLA_Obj_is_complex( A ) )
  {
    *gflops = ( 4.0 * (       1.0     * m * m * m ) + 
                4.0 * ( 4.0 / 3.0     * m * m * m ) + 
                      (       4.5 * k * m * m     ) +
                2.0 * (       3.0 * k * m * m * m ) ) / 
              dtime_old / 1e9;
  }
  else
  {
    *gflops = ( 1.0 * (       1.0     * m * m * m ) + 
                1.0 * ( 4.0 / 3.0     * m * m * m ) + 
                      (       4.5 * k * m * m     ) +
                1.0 * (       3.0 * k * m * m * m ) ) / 
              dtime_old / 1e9;
  }

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &Z );
}

