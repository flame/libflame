
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Symm_rl_unb_var9( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

  FLA_Obj CT,              C0,
          CB,              c1t,
                           C2;

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ){

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &b1t, 
                           BB,                &B2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &c1t, 
                           CB,                &C2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /* c1t  = c1t  + b1t * A    */
    /* c1t' = c1t' + A'  * b1t' */
    FLA_Symv_external( FLA_LOWER_TRIANGULAR, alpha, A, b1t, beta, c1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  b1t, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  c1t, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

#endif
