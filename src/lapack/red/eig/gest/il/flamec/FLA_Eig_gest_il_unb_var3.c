
#include "FLAME.h"

FLA_Error FLA_Eig_gest_il_unb_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
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

    // a10t = a10t - 1/2 * y10t;
    FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );

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

    // a10t = a10t - 1/2 * y10t;
    FLA_Axpy_external( FLA_MINUS_ONE_HALF, y10t, a10t );

    // a10t = inv(beta11) * a10t;
    FLA_Inv_scal_external( beta11, a10t );

    // Y20 = Y20 + b21 * a10t;
    FLA_Ger_external( FLA_ONE, b21, a10t, Y20 );

    // y21 = b21 * alpha11;
    FLA_Copy_external( b21, y21 );
    FLA_Scal_external( alpha11, y21 );

    // y21 = y21 + B20 * a10t';
    FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE,
                        FLA_ONE, B20, a10t, FLA_ONE, y21 );

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

