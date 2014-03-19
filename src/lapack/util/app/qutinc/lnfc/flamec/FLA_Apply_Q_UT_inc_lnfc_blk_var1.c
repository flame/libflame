/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_inc_lnfc_blk_var1( FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, fla_apqutinc_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TTL,   WTR,      T00, W01, W02, 
          TBL,   TBR,      T10, T11, W12,
                           T20, T21, T22;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t   b;
  dim_t   m_BR, n_BR;

  // If m > n, then we have to initialize our partitionings carefully so
  // that we begin in the proper location in A and B (since we traverse
  // matrix A from BR to TL).
  if ( FLA_Obj_length( A ) > FLA_Obj_width( A ) )
  {
    m_BR = FLA_Obj_length( A ) - FLA_Obj_width( A );
    n_BR = 0;
  }
  else
  {
    m_BR = 0;
    n_BR = 0;
  }

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     m_BR, n_BR, FLA_BR );

  FLA_Part_2x2( TW,   &TTL, &WTR,
                      &TBL, &TBR,     m_BR, n_BR, FLA_BR );

  FLA_Part_2x1( B,    &BT,
                      &BB,            m_BR, FLA_BOTTOM );

  while ( FLA_Obj_min_dim( ATL ) > 0 ){

    b = FLA_Determine_blocksize( ATL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x2_to_3x3( TTL, /**/ WTR,       &T00, &W01, /**/ &W02,
                                                &T10, &T11, /**/ &W12,
                        /* ************* */   /* ******************** */
                           TBL, /**/ TBR,       &T20, &T21, /**/ &T22,
                           b, b, FLA_TL );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                                              &B1,
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /*
        Apply Q^H to B1 and B2 from the left:
      
            / B1 \ = Q^H * / B1 \
            \ B2 /         \ B2 / 
      
        where Q is formed from A21 and T21.
    */

    FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              A21, T21, W1, B1,
                                            B2,
                              FLA_Cntl_sub_apq2ut( cntl ) );

    /*
        Apply Q^H to B1 from the left:
      
          B1 = Q^H * B1
      
        where Q is formed from A11 and T11. Note that W1 refers
        to a row-panel of blocks where each block refers to an
        nb_alg-by-b row-panel of workspace.
    */

    FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                             A11, T11, W1, B1,
                             FLA_Cntl_sub_apqut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &WTR,       T00, /**/ W01, W02,
                            /* ************** */  /* ****************** */
                                                     T10, /**/ T11, W12,
                              &TBL, /**/ &TBR,       T20, /**/ T21, T22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                            /* ** */           /* ** */
                                                  B1,
                              &BB,                B2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

