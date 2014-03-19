/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Eig_gest_il_unb_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  FLA_Obj yL,    yR,       y10t, psi11,  y12t;

  FLA_Obj y10t_t,
          y10t_b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_1x2( Y,    &yL,  &yR,      0, FLA_LEFT );

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

    FLA_Repart_1x2_to_1x3( yL,  /**/ yR,        &y10t, /**/ &psi11,  &y12t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( y10t,   &y10t_t,
                          &y10t_b,    1, FLA_TOP );

    // y10t   = 1/2 * b10t * A00;
    // y10t^T = 1/2 * A00^T * b10t^T;
    //        = 1/2 * conj(A00) * b10t^T;
    FLA_Hemvc_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
                        FLA_ONE_HALF, A00, b10t, FLA_ZERO, y10t_t );

    // a10t = a10t - y10t;
    FLA_Axpy_external( FLA_MINUS_ONE, y10t_t, a10t );

    // alpha11 = alpha11 - a10t * b10t' - b10t * a10t';
    FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a10t, b10t, FLA_ONE, alpha11 );

    // alpha11 = inv(beta11) * alpha11 * inv(conj(beta11));
    //         = inv(beta11) * alpha11 * inv(beta11);
    FLA_Inv_scal_external( beta11, alpha11 );
    FLA_Inv_scal_external( beta11, alpha11 );

    // a21 = a21 - A20 * b10t';
    FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
                        FLA_MINUS_ONE, A20, b10t, FLA_ONE, a21 );

    // a21 = a21 * inv(conj(beta11));
    //     = a21 * inv(beta11);
    FLA_Inv_scal_external( beta11, a21 );

    // a10t = a10t - y10t;
    FLA_Axpy_external( FLA_MINUS_ONE, y10t_t, a10t );

    // a10t = inv(beta11) * a10t;
    FLA_Inv_scal_external( beta11, a10t );

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

    FLA_Cont_with_1x3_to_1x2( &yL,  /**/ &yR,        y10t, psi11, /**/  y12t,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

