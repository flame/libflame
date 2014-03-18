
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Syr2k_un_unb_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

  FLA_Obj CTL,   CTR,      C00,  c01,     C02, 
          CBL,   CBR,      c10t, gamma11, c12t,
                           C20,  c21,     C22;


  FLA_Scalr_external( FLA_UPPER_TRIANGULAR, beta, C );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_BR );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){


    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &a1t, 
                        /* ** */            /* ** */
                           AB,                &A2,        1, FLA_TOP );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &b1t, 
                        /* ** */            /* ** */
                           BB,                &B2,        1, FLA_TOP );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00,  &c01,     /**/ &C02,
                                                &c10t, &gamma11, /**/ &c12t,
                        /* ************* */   /* ************************** */
                           CBL, /**/ CBR,       &C20,  &c21,     /**/ &C22,
                           1, 1, FLA_TL );

    /*------------------------------------------------------------*/

    /* c01 = c01 + A0 * b1t' */
    FLA_Gemv_external( FLA_NO_TRANSPOSE, alpha, A0, b1t, FLA_ONE, c01 );    

    /* c12t = c12t + b1t * A2' */
    FLA_Gemv_external( FLA_NO_TRANSPOSE, alpha, A2, b1t, FLA_ONE, c12t );

    /* gamma11 = gamma11 + a1t * b1t' + b1t * a1t' */
    FLA_Dot2s_external( alpha, a1t, b1t, FLA_ONE, gamma11 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  a1t, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  b1t, 
                              &BB,                B2,     FLA_BOTTOM );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00,  /**/ c01,     C02,
                            /* ************** */  /* ************************ */
                                                     c10t, /**/ gamma11, c12t,
                              &CBL, /**/ &CBR,       C20,  /**/ c21,     C22,
                              FLA_BR );

  }

  return FLA_SUCCESS;
}

#endif
