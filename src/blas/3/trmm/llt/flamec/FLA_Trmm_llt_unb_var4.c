
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trmm_llt_unb_var4( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Scal_external( alpha, B );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( BR ) < FLA_Obj_width( B ) ){

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &b1, /**/ &B2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* b1 = tril( A )' * b1 */
    FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diagA, A, b1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ b1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
