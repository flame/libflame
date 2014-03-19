/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_LU_piv_blk_var4( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  FLA_Error r_val = FLA_SUCCESS, r_val_sub = FLA_SUCCESS;
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj AB0, AB1, AB2;

  dim_t b;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
          FLA_Obj_width( ATL ) < FLA_Obj_width( A )){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &p1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // A11 = A11 - A10 * A0
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A01, FLA_ONE, A11,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // A21 = A21 - A20 * A01
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A20, A01, FLA_ONE, A21,
                       FLA_Cntl_sub_gemm3( cntl ) );

    // AB1 = / A11 \
    //       \ A21 /
    FLA_Merge_2x1( A11,
                   A21,      &AB1 );

    // AB1, p1 = LU_piv( AB1 )
    FLA_LU_piv_internal( AB1, p1, 
                         FLA_Cntl_sub_lu( cntl ) );

    // If the unblocked algorithm returns a null pivot, 
    // update the pivot index and return it.
    if ( r_val == FLA_SUCCESS && r_val_sub >= 0 )
    {
        r_val = FLA_Obj_length( A01 ) + r_val_sub;
    }

    // AB0 = / A10 \
    //       \ A20 /
    FLA_Merge_2x1( A10,
                   A20,      &AB0 );

    // AB2 = / A12 \
    //       \ A22 /
    FLA_Merge_2x1( A12,
                   A22,      &AB2 );

    // Apply pivots to remaining columns
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB0,
                               FLA_Cntl_sub_appiv1( cntl ) );
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB2,
                               FLA_Cntl_sub_appiv1( cntl ) );

    // A12 = A12 - A10 * A02
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A02, FLA_ONE, A12,
                       FLA_Cntl_sub_gemm2( cntl ) );

    // A12 = trilu( A11 ) \ A12 
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, A12,
                       FLA_Cntl_sub_trsm1( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  p1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );

  }

  return r_val;
}

#endif
