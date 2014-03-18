
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Hemm_rl_blk_var10( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl )
{
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  dim_t b;

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( BB ) < FLA_Obj_length( B ) ){

    b = FLA_Determine_blocksize( BT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &C1, 
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* C1 = C1 + B1 * A */
    FLA_Hemm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR, 
                       alpha, A, B1, beta, C1,
                       FLA_Cntl_sub_hemm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* ** */
                                                  C1, 
                              &CB,                C2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
