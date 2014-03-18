
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Gemm_tn_unb_var4( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Obj CL,    CR,       C0,  c1,  C2;

  FLA_Scal_external( beta, C );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( BR ) < FLA_Obj_width( B ) ){

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &b1, /**/ &B2,
                           1, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, &c1, /**/ &C2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* c1 = A' * b1 + c1 */
    FLA_Gemv_external( FLA_TRANSPOSE, alpha, A, b1, FLA_ONE, c1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ b1, B2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, /**/ c1, C2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
