
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_tn_unb_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj CT,              C0,
          CB,              c1t,
                           C2;

  FLA_Scal_external( beta, C );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &a1, /**/ &A2,
                           1, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &c1t, 
                        /* ** */            /* *** */
                           CB,                &C2,        1, FLA_TOP );

    /*------------------------------------------------------------*/

    /* c1t  = a1' * B + c1t  */
    /* c1t' = B' * a1 + c1t' */
    FLA_Gemv_external( FLA_TRANSPOSE, alpha, B, a1, FLA_ONE, c1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ a1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* *** */
                                                  c1t, 
                              &CB,                C2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
