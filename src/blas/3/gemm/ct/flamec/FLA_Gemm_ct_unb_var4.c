
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_ct_unb_var4( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

  FLA_Obj CL,    CR,       C0,  c1,  C2;

  FLA_Scal_external( beta, C );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( BB ) < FLA_Obj_length( B ) ){

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &b1t, 
                        /* ** */            /* *** */
                           BB,                &B2,        1, FLA_TOP );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, &c1, /**/ &C2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* c1 = A * b1t + c1 */
    FLA_Gemv_external( FLA_CONJ_NO_TRANSPOSE, alpha, A, b1t, FLA_ONE, c1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* *** */
                                                  b1t, 
                              &BB,                B2,     FLA_BOTTOM );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, /**/ c1, C2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
