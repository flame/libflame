
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trsm_ruh_unb_var4( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

  FLA_Scal_external( alpha, B );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( BB ) < FLA_Obj_length( B ) ){

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &b1t, 
                        /* ** */            /* *** */
                           BB,                &B2,        1, FLA_TOP );

    /*------------------------------------------------------------*/

    /* b1t = b1t / triu( A' ); */
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diagA, A, b1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* *** */
                                                  b1t, 
                              &BB,                B2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
