
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Herk_un_unb_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj CTL,   CTR,      C00,  c01,     C02, 
          CBL,   CBR,      c10t, gamma11, c12t,
                           C20,  c21,     C22;

  FLA_Scalr_external( FLA_UPPER_TRIANGULAR, beta, C );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00,  /**/ &c01,     &C02,
                        /* ************* */   /* ************************** */
                                                &c10t, /**/ &gamma11, &c12t,
                           CBL, /**/ CBR,       &C20,  /**/ &c21,     &C22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    /* c01     = c01     + A0  * a1t' */
    FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE, alpha, A0, a1t, FLA_ONE, c01 );

    /* gamma11 = gamma11 + a1t * a1t' */
    FLA_Dotcs_external( FLA_CONJUGATE, alpha, a1t, a1t, FLA_ONE, gamma11 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00,  c01,     /**/ C02,
                                                     c10t, gamma11, /**/ c12t,
                            /* ************** */  /* ************************ */
                              &CBL, /**/ &CBR,       C20,  c21,     /**/ C22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}

#endif
