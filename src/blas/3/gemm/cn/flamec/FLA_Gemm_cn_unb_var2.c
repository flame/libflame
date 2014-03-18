
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_cn_unb_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj CT,              C0,
          CB,              c1t,
                           C2;

  FLA_Scal_external( beta, C );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &a1t, 
                        /* ** */            /* *** */
                           AB,                &A2,        1, FLA_TOP );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &c1t, 
                        /* ** */            /* *** */
                           CB,                &C2,        1, FLA_TOP );

    /*------------------------------------------------------------*/

    /* c1t  = a1t * B + c1t    */
    /* c1t' = B' * a1t' + c1t' */
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, alpha, B, a1t, FLA_ONE, c1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* *** */
                                                  a1t, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* *** */
                                                  c1t, 
                              &CB,                C2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
