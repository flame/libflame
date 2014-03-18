
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_ht_unb_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Scal_external( beta, C );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &a1t, 
                        /* ** */            /* *** */
                           AB,                &A2,        1, FLA_TOP );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &b1, /**/ &B2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* C = a1t' * b1' + C */
    FLA_Gerc_external( FLA_CONJUGATE, FLA_NO_CONJUGATE, alpha, a1t, b1, C );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* *** */
                                                  a1t, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ b1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
