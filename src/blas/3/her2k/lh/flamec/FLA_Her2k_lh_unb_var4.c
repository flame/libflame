
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Her2k_lh_unb_var4( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Obj CTL,   CTR,      C00,  c01,     C02, 
          CBL,   CBR,      c10t, gamma11, c12t,
                           C20,  c21,     C22;


  FLA_Scalr_external( FLA_LOWER_TRIANGULAR, beta, C );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){


    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &b1, &B2,
                           1, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00,  /**/ &c01,     &C02,
                        /* ************* */   /* ************************** */
                                                &c10t, /**/ &gamma11, &c12t,
                           CBL, /**/ CBR,       &C20,  /**/ &c21,     &C22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    /* c21 = c21 + A2' * b1 */
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, alpha, A2, b1, FLA_ONE, c21 );    

    /* c21 = c21 + B2' * a1 */
    FLA_Gemv_external( FLA_CONJ_TRANSPOSE, alpha, B2, a1, FLA_ONE, c21 );    

    /* gamma11 = gamma11 + a1' * b1 + b1' * a1 */
    FLA_Dot2cs_external( FLA_CONJUGATE, alpha, a1, b1, FLA_ONE, gamma11 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, b1, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00,  c01,     /**/ C02,
                                                     c10t, gamma11, /**/ c12t,
                            /* ************** */  /* ************************ */
                              &CBL, /**/ &CBR,       C20,  c21,     /**/ C22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

#endif
