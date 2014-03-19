/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

// w should already has colnorms.
FLA_Error FLA_QR_UT_piv_unb_var1( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p )
{
  FLA_Obj AL,    AR;

  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;
  FLA_Obj AB1, AT1, at1;

  FLA_Obj TTL,   TTR,      T00,  t01,   T02, 
          TBL,   TBR,      t10t, tau11, t12t,
                           T20,  t21,   T22;

  FLA_Obj pT,              p0,
          pB,              pi1,
                           p2;

  FLA_Obj wT,              w0,
          wB,              omega1,
                           w2;

  dim_t   nb  = FLA_Obj_width ( A ) - FLA_Obj_width( T );
  //dim_t   mb  = FLA_Obj_length( A ) - FLA_Obj_width( T );

  FLA_Part_1x2( A,    &AL,  &AR,      nb, FLA_RIGHT );

  FLA_Part_2x2( AL,   &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( AR,   &AT,
                      &AB,            0, FLA_TOP );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  FLA_Part_2x1( p,    &pT,
                      &pB,            0, FLA_TOP );

  FLA_Part_2x1( w,    &wT,
                      &wB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* *** */
                                              &a1t,
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( pT,                &p0,
                        /* ** */            /* *** */
                                              &pi1,
                           pB,                &p2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( wT,                &w0,
                        /* ** */            /* *** */
                                              &omega1,
                           wB,                &w2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/


    // Ignore minus inputs for LAPACK compatability.
    if ( FLA_Obj_le( pi1, FLA_ZERO ) == FALSE )
    {
      // Determine pivot index
      FLA_Amax_external( wB, pi1 );
      FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, wB );
      
      // Apply pivots
      FLA_Merge_1x2( ABR, AB, &AB1 );
      FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, pi1, AB1 );
    }
    else
    {
      // Do not pivot.
      FLA_Set( FLA_ZERO, pi1 );
    }

    // Compute tau11 and u21 from alpha11 and a21 such that tau11 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    FLA_Househ2_UT( FLA_LEFT,
                    alpha11,
                    a21, tau11 );

    // Apply H to (a12t A22)^T
    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau11 and u21.
    FLA_Merge_1x2( A22,  A2,  &AB1 );
    FLA_Merge_1x2( a12t, a1t, &at1 );

    FLA_Apply_H2_UT( FLA_LEFT, tau11, a21, at1,
                                           AB1 );

    // t01 = a10t' + A20' * u21;
    FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, t01 );
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, t01 );

    // Apply pivots to previous rows
    if ( FLA_Obj_le( pi1, FLA_ZERO ) == FALSE )
    {
      FLA_Merge_1x2( ATR, AT, &AT1 );
      FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, pi1, AT1 );
    }

    // Norm downdate w2 = alpha w2 + beta columnwisenorm2(a12t)
    FLA_QR_UT_piv_colnorm( FLA_MINUS_ONE, at1, w2 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  a1t,
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0,
                                                  pi1,
                            /* ** */           /* *** */
                              &pB,                p2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &wT,                w0,
                                                  omega1,
                            /* ** */           /* *** */
                              &wB,                w2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

