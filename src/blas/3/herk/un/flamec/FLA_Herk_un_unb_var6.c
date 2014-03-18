
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Herk_un_unb_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Scalr_external( FLA_UPPER_TRIANGULAR, beta, C );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &a1, /**/ &A2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* C := C + a1 * a1' */ 
    FLA_Her_external( FLA_UPPER_TRIANGULAR, alpha, a1, C );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ a1, A2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
