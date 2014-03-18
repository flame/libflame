
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Hemm_rl_unb_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BL,    BR,       B0,  b1t,  B2;

  FLA_Obj CL,    CR,       C0,  c1t,  C2;

  FLA_Scal_external( beta, C );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &b1t, &B2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &c1t, &C2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* c1t = c1t + B0 * a10t' */
    FLA_Gemvc_external( FLA_NO_TRANSPOSE, FLA_CONJUGATE, alpha, B0, a10t, FLA_ONE, c1t );

    /* c1t = c1t + b1t * alpha11 */
    FLA_Axpys_external( alpha, alpha11, b1t, FLA_ONE, c1t );

    /* c1t = c1t + B2 * a21 */
    FLA_Gemv_external( FLA_NO_TRANSPOSE, alpha, B2, a21, FLA_ONE, c1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, b1t, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, c1t, /**/ C2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

#endif
