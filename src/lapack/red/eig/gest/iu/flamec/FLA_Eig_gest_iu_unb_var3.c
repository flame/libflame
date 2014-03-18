
#include "FLAME.h"

FLA_Error FLA_Eig_gest_iu_unb_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  FLA_Obj YTL,   YTR,      Y00,  y01,   Y02,
          YBL,   YBR,      y10t, psi11, y12t,
                           Y20,  y21,   Y22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x2( Y,    &YTL, &YTR,
                      &YBL, &YBR,     0, 0, FLA_TL );

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

    FLA_Repart_2x2_to_3x3( YTL, /**/ YTR,       &Y00,  /**/ &y01,   &Y02,
                        /* ************* */   /* ************************ */
                                                &y10t, /**/ &psi11, &y12t,
                           YBL, /**/ YBR,       &Y20,  /**/ &y21,   &Y22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    // a01 = a01 - 1/2 * y01;
    FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );

    // alpha11 = alpha11 - a01' * b01 - b01' * a01;
    FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );

    // alpha11 = inv(beta11) * alpha11 * inv(conj(beta11));
    //         = inv(beta11) * alpha11 * inv(beta11);
    FLA_Inv_scal_external( beta11, alpha11 );
    FLA_Inv_scal_external( beta11, alpha11 );

    // a12t   = a12t   - b01' * A02; 
    // a12t^T = a12t^T - A02^T * conj(b01); 
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
                        FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );

    // a12t = inv(conj(beta11)) * a12t;
    // a12t = inv(beta11) * a12t;
    FLA_Inv_scal_external( beta11, a12t );

    // a01 = a01 - 1/2 * y01;
    FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );

    // a01 = a01 * inv(beta11);
    FLA_Inv_scal_external( beta11, a01 );

    // Y02 = Y02 + a01 * b12t;
    FLA_Ger_external( FLA_ONE, a01, b12t, Y02 );

    // y12t = alpha11 * b12t;
    FLA_Copy_external( b12t, y12t );
    FLA_Scal_external( alpha11, y12t );

    // y12t   = y12t   + a01' * B02;
    // y12t^T = y12t^T + B02^T * conj(a01);
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
                        FLA_ONE, B02, a01, FLA_ONE, y12t );

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

    FLA_Cont_with_3x3_to_2x2( &YTL, /**/ &YTR,       Y00,  y01,   /**/ Y02,
                                                     y10t, psi11, /**/ y12t,
                            /* ************** */  /* ********************** */
                              &YBL, /**/ &YBR,       Y20,  y21,   /**/ Y22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

