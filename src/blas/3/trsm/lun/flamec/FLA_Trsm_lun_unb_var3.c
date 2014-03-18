
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trsm_lun_unb_var3( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Scal_external( alpha, B );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &b1, &B2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* b1 = triu( A ) \ b1; */
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diagA, A, b1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, b1, /**/ B2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

#endif
