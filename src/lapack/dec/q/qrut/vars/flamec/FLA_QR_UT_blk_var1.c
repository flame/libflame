/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR_UT_blk_var1( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  W12;

  FLA_Obj T1T,   T2B;

  FLA_Obj AB1,   AB2;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = fla_min( b_alg, FLA_Obj_min_dim( ABR ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,   &T1T, 
                        &T2B,    b, FLA_TOP );

    FLA_Merge_2x1( A11,
                   A21,   &AB1 );

    // Perform a QR factorization via the UT transform on AB1:
    //
    //   / A11 \ -> QB1 R11
    //   \ A21 /
    //
    // where:
    //  - QB1 is formed from UB1 (which is stored column-wise below the
    //    diagonal of AB1) and T11 (which is stored to the upper triangle
    //    of T11).
    //  - R11 is stored to the upper triangle of AB1.
  
    FLA_QR_UT_internal( AB1, T1T, 
                        FLA_Cntl_sub_qrut( cntl ) );


    if ( FLA_Obj_width( A12 ) > 0 )
    {
      FLA_Merge_2x1( A12,
                     A22,   &AB2 );

      // Apply the Householder transforms associated with UB1 and T11 to 
      // AB2:
      //
      //   / A12 \ := QB1' / A12 \
      //   \ A22 /         \ A22 /
      //
      // where QB1 is formed from UB1 and T11.

      FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                               AB1, T1T, W12, AB2,
                               FLA_Cntl_sub_apqut( cntl ) );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

