/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Eig_gest_nu_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  //FLA_Obj yT,              y01,
  //        yB,              psi11,
  //                         y21;

  //FLA_Obj y01_l, y01_r;

  FLA_Obj psi11, y12t,
          y21,   Y22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  //FLA_Part_2x1( Y,    &yT,
  //                    &yB,            0, FLA_TOP );

  FLA_Part_2x2( Y,    &psi11, &y12t,
                      &y21,   &Y22,     1, 1, FLA_TL );

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

    //FLA_Repart_2x1_to_3x1( yT,                  &y01,
    //                    /* ** */              /* ***** */
    //                                            &psi11,
    //                       yB,                  &y21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    //FLA_Part_1x2( y01,    &y01_l, &y01_r,     1, FLA_LEFT );

    //// y01 = b01 * alpha11;
    //FLA_Copy_external( b01, y01_l );
    //FLA_Scal_external( alpha11, y01_l );
    // psi11 = 1/2 * alpha11;
    FLA_Copy_external( alpha11, psi11 );
    FLA_Scal_external( FLA_ONE_HALF, psi11 );

    // a01 = triu( B00 ) * a01;
    FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       B00, a01 );

    //// a01 = a01 + 1/2 * y01;
    //FLA_Axpy_external( FLA_ONE_HALF, y01_l, a01 );
    // a01 = a01 + 1/2 * alpha11 * y01;
    FLA_Axpy_external( psi11, b01, a01 );

    // A00 = A00 + a01 * b01' + b01 * a01';
    FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
                        FLA_ONE, a01, b01, A00 );

    //// a01 = a01 + 1/2 * y01;
    //FLA_Axpy_external( FLA_ONE_HALF, y01_l, a01 );
    // a01 = a01 + 1/2 * alpha11 * y01;
    FLA_Axpy_external( psi11, b01, a01 );

    // a01 = a01 * beta11;
    FLA_Scal_external( beta11, a01 );

    // alpha11 = conj(beta11) * alpha11 * beta11;
    //         = beta11 * alpha11 * beta11;
    FLA_Scal_external( beta11, alpha11 );
    FLA_Scal_external( beta11, alpha11 );

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

    //FLA_Cont_with_3x1_to_2x1( &yT,                   y01,
    //                                                 psi11,
    //                        /* ** */              /* ***** */
    //                          &yB,                   y21,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

#endif
