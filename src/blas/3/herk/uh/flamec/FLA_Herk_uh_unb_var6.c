
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Herk_uh_unb_var6( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Scalr_external( FLA_UPPER_TRIANGULAR, beta, C );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &a1t, 
                        /* ** */            /* *** */
                           AB,                &A2,        1, FLA_TOP );

    /*------------------------------------------------------------*/

    /* C := C + a1t' * a1t */ 
    FLA_Herc_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE, alpha, a1t, C );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* *** */
                                                  a1t, 
                              &AB,                A2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
