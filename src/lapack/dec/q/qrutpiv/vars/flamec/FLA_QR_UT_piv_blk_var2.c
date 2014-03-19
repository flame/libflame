/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_blk_var2( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  W12;
  FLA_Obj TT,    TB;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj wT,              w0,
          wB,              w1,
                           w2;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  FLA_Part_2x1( w,    &wT, 
                      &wB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( ABR ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &p1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( wT,                &w0, 
                        /* ** */            /* ** */
                                              &w1, 
                           wB,                &w2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // ** Reshape T matrices to match the blocksize b
    FLA_Part_2x1( TR,   &TT, 
                        &TB,    b, FLA_TOP );

    // ** Perform a unblocked (BLAS2-oriented) QR factorization 
    // with pivoting via the UT transform on ABR:
    //
    //   ABR  -> QB1 R11
    //
    // where:
    //  - QB1 is formed from UB1 (which is stored column-wise below the
    //    diagonal of ( A11 A21 )^T and the upper-triangle of T1. 
    //  - R11 is stored to ( A11 A12 ).
    //  - W12 stores  T and partial updates for FLA_Apply_Q_UT_piv_var.
    FLA_QR_UT_piv_internal( ABR, TT, wB, p1, 
                            FLA_Cntl_sub_qrut( cntl ) );

    if ( FLA_Obj_width( A12 ) > 0 )
    {
      // ** Block update
      FLA_Part_2x1( W12,  &TT, 
                          &TB,    b, FLA_TOP );
 
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                         FLA_MINUS_ONE, A21, TT, FLA_ONE, A22 );
    }

    // ** Apply pivots to previous columns.
    FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, p1, ATR );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  p1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &wT,                w0, 
                                                  w1, 
                            /* ** */           /* ** */
                              &wB,                w2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

