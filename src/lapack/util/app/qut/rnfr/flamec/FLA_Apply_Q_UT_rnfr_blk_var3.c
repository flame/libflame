/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_rnfr_blk_var3( FLA_Obj A, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TWTL,  TWTR,     TW00, TW01, TW02, 
          TWBL,  TWBR,     TW10,  T11,  W12,
                           TW20, TW21, TW22;

  FLA_Obj WTL,   WTR,
          WBL,   WBR;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( TW,   &TWTL, &TWTR,
                      &TWBL, &TWBR,     0, 0, FLA_TL );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TWTL, /**/ TWTR,       &TW00, /**/ &TW01, &TW02,
                        /* *************** */   /* *********************** */
                                                  &TW10, /**/  &T11,  &W12,
                           TWBL, /**/ TWBR,       &TW20, /**/ &TW21, &TW22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( W,    &WTL, &WTR,
                        &WBL, &WBR,     b, FLA_Obj_length( B1 ), FLA_TL );

    // WTL = B1;

    FLA_Copyt_internal( FLA_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = trilu( A11 );
    // U12 = A12;
    // Let WTL^T be conformal to B1.
    //
    // WTL^T = ( B1 * U11^T + B2 * U12^T ) * inv( triu(T11) );
    // WTL   = inv( triu(T11) )^T * ( U11 * B1^T + U12 * B2^T );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, 
                       FLA_ONE, A12, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T11, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - WTL^T * conj(U12);
    // B1 = B1 - WTL^T * conj(U11);
    //    = B1 - ( U11' * WTL )^T;

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE,
                       FLA_MINUS_ONE, WTL, A12, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TWTL, /**/ &TWTR,       TW00, TW01, /**/ TW02,
                                                       TW10,  T11, /**/  W12,
                            /* **************** */  /* ********************* */
                              &TWBL, /**/ &TWBR,       TW20, TW21, /**/ TW22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

