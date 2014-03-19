/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT1  3

/*
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

extern fla_apqut_t* fla_apqut_cntl_leaf;
*/

FLA_Error REF_Apply_Q_UT_lnfc( FLA_Obj A, FLA_Obj t, FLA_Obj B );
FLA_Error REF_Bidiag_form_U_blk_external( FLA_Side side, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );
void time_Apply_Q_UT_lnfc(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj A_orig, FLA_Obj t, FLA_Obj T, FLA_Obj s, FLA_Obj S, FLA_Obj B,
               double *dtime, double *diff, double *gflops );


void time_Apply_Q_UT_lnfc(
               int variant, int type, int n_repeats, int m, int n, int nb_alg,
               FLA_Obj A, FLA_Obj A_orig, FLA_Obj t, FLA_Obj T, FLA_Obj s, FLA_Obj S, FLA_Obj B,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, A_orig_save, B_save, norm;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_orig_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &B_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( A, A_orig_save );
  FLA_Copy_external( B, B_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );
    FLA_Copy_external( A_orig_save, A_orig );
    FLA_Copy_external( B_save, B );

    *dtime = FLA_Clock();

    switch( variant )
    {

    case 0:
      REF_Apply_Q_UT_lnfc( A, t, B );
      //REF_Bidiag_form_U_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, A, t, B );
      //FLA_Bidiag_blk_external( A_orig, t, s );
      //REF_Bidiag_form_U_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, A_orig, t, B );
      break;

    case 1:
    {
      // Time variant 1 
      switch( type ){
      case FLA_ALG_BLOCKED:
        //FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, B );
        FLA_QR_UT_form_Q( A, T, B );
        //FLA_Bidiag_UT_form_U( A, T, B );
        //FLA_Bidiag_UT( A_orig, T, S );
        //FLA_Bidiag_UT_form_U( A_orig, T, B );
        break;
      }

      break;
    }


    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }



/*
  if ( variant == 0 )
  {
    FLA_Copy_external( b, b_ref );
    if ( FLA_Obj_is_real( A ) )
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_TRANSPOSE, FLA_COLUMNWISE, A, t, b );
    else
      FLA_Apply_Q_blk_external( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_COLUMNWISE, A, t, b );
    FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    if ( FLA_Obj_is_single_precision( A ) )
      *diff = *(FLA_FLOAT_PTR(norm));
    else
      *diff = *(FLA_DOUBLE_PTR(norm));
  }
  else
*/
  {
    FLA_Obj_set_to_identity( A );
//FLA_Obj_show( "B", B, "%8.1e %8.1e ", "" );
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, B, B, FLA_MINUS_ONE, A );
    FLA_Norm_frob( A, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

  }
/*
  *gflops = 2.0 * n * n * ( m - n / 3.0 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;
*/

    *gflops = ( 4.0 * ( 2.0           * m * n * n - 2.0 / 3.0 * n * n * n ) +
                4.0 * ( 4.0 / 3.0     * m * m * m ) +
                4.0 * ( 4.0 / 3.0     * n * n * n ) +
                      ( 13.0      * 2 * m * m     ) +
                2.0 * (       3.0 * 2 * m * m * m ) +
                2.0 * (       3.0 * 2 * n * n * n ) ) /
              dtime_old / 1e9;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( A_orig_save, A_orig );
  FLA_Copy_external( B_save, B );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &A_orig_save );
  FLA_Obj_free( &B_save );
  FLA_Obj_free( &norm );
}

