/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Lyap_n_blk_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;
  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;
  dim_t   b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_BR );

  while ( FLA_Obj_length( CTL ) > 0 ){

    b = FLA_Determine_blocksize( CTL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, &C01, /**/ &C02,
                                                &C10, &C11, /**/ &C12,
                        /* ************* */   /* ******************** */
                           CBL, /**/ CBR,       &C20, &C21, /**/ &C22,
                           b, b, FLA_TL );

    /*------------------------------------------------------------*/

    // C12 = sylv( A11, A22', C12 );
    FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, A11, A22, C12, scale,
                       FLA_Cntl_sub_sylv( cntl ) );
  
    // C11 = isgn * C11 - A12 * C12' - C12 * A12';
    // C11 = lyap_n( A11, C11 );
    FLA_Her2k_internal( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                        FLA_MINUS_ONE, A12, C12, isgn, C11,
                        FLA_Cntl_sub_her2k( cntl ) );
    FLA_Lyap_internal( FLA_NO_TRANSPOSE, FLA_ONE, A11, C11, scale,
                       FLA_Cntl_sub_lyap( cntl ) );

    // C02 = C02 - A01 * C12;
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A01, C12, FLA_ONE, C02,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // C01 = isgn * C01 - A01 * C11;
    // C01 = C01 - A02 * C12';
    FLA_Hemm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                       FLA_MINUS_ONE, C11, A01, isgn, C01,
                       FLA_Cntl_sub_hemm( cntl ) );
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_MINUS_ONE, A02, C12, FLA_ONE, C01,
                       FLA_Cntl_sub_gemm2( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, /**/ C01, C02,
                            /* ************** */  /* ****************** */
                                                     C10, /**/ C11, C12,
                              &CBL, /**/ &CBR,       C20, /**/ C21, C22,
                              FLA_BR );
  }

  return FLA_SUCCESS;
}

