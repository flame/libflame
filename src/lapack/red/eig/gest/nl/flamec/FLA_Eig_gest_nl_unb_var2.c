/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Eig_gest_nl_unb_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  FLA_Obj yT,              y01,
          yB,              psi11,
                           y21;

  FLA_Obj y21_l, y21_r;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x1( Y,    &yT,
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00,  /**/ &b01,    &B02,
                        /* ************* */   /* ************************* */
                                                &b10t, /**/ &beta11, &b12t,
                           BBL, /**/ BBR,       &B20,  /**/ &b21,    &B22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( yT,                  &y01,
                        /* ** */              /* ***** */
                                                &psi11,
                           yB,                  &y21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( y21,    &y21_l, &y21_r,     1, FLA_LEFT );

    // a10t = beta11 * a10t;
    FLA_Scal_external( beta11, a10t );

    // a10t   = a10t   + b21' * A20;
    // a10t^T = a10t^T + A20^T * conj(b21);
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
                        FLA_ONE, A20, b21, FLA_ONE, a10t );

    // y21 = A22 * b21;
    FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
                       FLA_ONE, A22, b21, FLA_ZERO, y21_l );

    // a21 = a21 * beta11;
    FLA_Scal_external( beta11, a21 );

    // a21 = a21 + 1/2 * y21;
    FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );

    // alpha11 = conj(beta11) * alpha11 * beta11;
    //         = beta11 * alpha11 * beta11;
    FLA_Scal_external( beta11, alpha11 );
    FLA_Scal_external( beta11, alpha11 );

    // alpha11 = alpha11 + a21' * b21 + b21' * a21;
    FLA_Dot2cs_external( FLA_CONJUGATE, FLA_ONE, a21, b21, FLA_ONE, alpha11 );

    // a21 = a21 + 1/2 * y21;
    FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00,  b01,    /**/ B02,
                                                     b10t, beta11, /**/ b12t,
                            /* ************** */  /* *********************** */
                              &BBL, /**/ &BBR,       B20,  b21,    /**/ B22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &yT,                   y01,
                                                     psi11,
                            /* ** */              /* ***** */
                              &yB,                   y21,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

