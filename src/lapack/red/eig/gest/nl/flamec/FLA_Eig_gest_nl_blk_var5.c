/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Eig_gest_nl_blk_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02,
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj BTL,   BTR,      B00, B01, B02,
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  FLA_Obj YTL,   YTR,      Y00, Y01, Y02,
          YBL,   YBR,      Y10, Y11, Y12,
                           Y20, Y21, Y22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x2( Y,    &YTL, &YTR,
                      &YBL, &YBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00, /**/ &B01, &B02,
                        /* ************* */   /* ******************** */
                                                &B10, /**/ &B11, &B12,
                           BBL, /**/ BBR,       &B20, /**/ &B21, &B22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( YTL, /**/ YTR,       &Y00, /**/ &Y01, &Y02,
                        /* ************* */   /* ******************** */
                                                &Y10, /**/ &Y11, &Y12,
                           YBL, /**/ YBR,       &Y20, /**/ &Y21, &Y22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    // Y10 = A11 * B10;
    FLA_Hemm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_ONE, A11, B10, FLA_ZERO, Y10,
                       FLA_Cntl_sub_hemm( cntl ) );

    // A10 = A10 * tril( B00 );
    FLA_Trmm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B00, A10,
                       FLA_Cntl_sub_trmm1( cntl ) );

    // A10 = A10 + 1/2 * Y10;
    FLA_Axpy_internal( FLA_ONE_HALF, Y10, A10,
                       FLA_Cntl_sub_axpy1( cntl ) );

    // A00 = A00 + A10' * B10 + B10' * A10;
    FLA_Her2k_internal( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                        FLA_ONE, A10, B10, FLA_ONE, A00,
                        FLA_Cntl_sub_her2k( cntl ) );

    // A10 = A10 + 1/2 * Y10;
    FLA_Axpy_internal( FLA_ONE_HALF, Y10, A10,
                       FLA_Cntl_sub_axpy2( cntl ) );

    // A10 = tril( B11 )' * A10;
    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B11, A10,
                       FLA_Cntl_sub_trmm2( cntl ) );

    // A11 = tril( B11 )' * A11 * tril( B11 );
    FLA_Eig_gest_internal( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR,
                           A11, Y11, B11,
                           FLA_Cntl_sub_eig_gest( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00, B01, /**/ B02,
                                                     B10, B11, /**/ B12,
                            /* ************** */  /* ****************** */
                              &BBL, /**/ &BBR,       B20, B21, /**/ B22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &YTL, /**/ &YTR,       Y00, Y01, /**/ Y02,
                                                     Y10, Y11, /**/ Y12,
                            /* ************** */  /* ****************** */
                              &YBL, /**/ &YBR,       Y20, Y21, /**/ Y22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

#endif
