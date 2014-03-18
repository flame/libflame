
#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_rnfr_blk_var2( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  FLA_Obj WL,    WR,       W0,  W1,  W2;

  dim_t b;

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  FLA_Part_1x2( W,    &WL,  &WR,      0, FLA_LEFT );

  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ){

    b = FLA_Determine_blocksize( BT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &B1, 
                           BB,                &B2,        b, FLA_BOTTOM );

    FLA_Repart_1x2_to_1x3( WL,  /**/ WR,        &W0, /**/ &W1, &W2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // B1 = B1 * Q;
    FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                             A, T, W1, B1,
                             FLA_Cntl_sub_apqut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  B1, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );

    FLA_Cont_with_1x3_to_1x2( &WL,  /**/ &WR,        W0, W1, /**/ W2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

