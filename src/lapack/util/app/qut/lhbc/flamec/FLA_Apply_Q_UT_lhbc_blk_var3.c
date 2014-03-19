/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_lhbc_blk_var3( FLA_Obj A, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TWTL,  TWTR,     TW00, TW01, TW02, 
          TWBL,  TWBR,     TW10,  T11,  W12,
                           TW20, TW21, TW22;

  FLA_Obj WTL,  WTR,
          WBL,  WBR;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t b, m_BR, n_BR;

  // If m < n, then we have to initialize our partitionings carefully so
  // that we begin in the proper location in A and B (since we traverse
  // matrix A from BR to TL).
  if ( FLA_Obj_length( A ) < FLA_Obj_width( A ) )
  {
    m_BR = 0;
    n_BR = FLA_Obj_width( A ) - FLA_Obj_length( A );
  }
  else
  {
    m_BR = 0;
    n_BR = 0;
  }

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     m_BR, n_BR, FLA_BR );

  FLA_Part_2x2( TW,   &TWTL, &TWTR,
                      &TWBL, &TWBR,    0, 0, FLA_BR );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            n_BR, FLA_BOTTOM );

  while ( FLA_Obj_min_dim( ATL ) > 0 ){

    b = FLA_Determine_blocksize( ATL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    // Since T was filled from TL to BR, and since we need to access them in
    // reverse order, we need to handle the case where the bottom-right-most
    // block is smaller than the other b x b blocks.
    if ( FLA_Obj_width( TWBR ) == 0 && 
         FLA_Obj_width( TW ) % b > 0 )
      b = FLA_Obj_width( TW ) % b;

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x2_to_3x3( TWTL, /**/ TWTR,       &TW00, &TW01, /**/ &TW02,
                                                  &TW10,  &T11, /**/  &W12,
                        /* *************** */   /* *********************** */
                           TWBL, /**/ TWBR,       &TW20, &TW21, /**/ &TW22,
                           b, b, FLA_TL );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_width( B1 ), FLA_TL );

    // WTL = B1;

    FLA_Copyt_internal( FLA_NO_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = trilu( A11 );
    // U21 = A21;
    // 
    // WTL = inv( triu(T11) ) * ( U11' * B1 + U21' * B2 );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, A21, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T11, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - U21 * WTL;
    // B1 = B1 - U11 * WTL;

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A21, WTL, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_NO_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x3_to_2x2( &TWTL, /**/ &TWTR,       TW00, /**/ TW01, TW02,
                            /* **************** */  /* ********************* */
                                                       TW10, /**/  T11,  W12,
                              &TWBL, /**/ &TWBR,       TW20, /**/ TW21, TW22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

