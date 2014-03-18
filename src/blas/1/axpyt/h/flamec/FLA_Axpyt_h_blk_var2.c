
#include "FLAME.h"

FLA_Error FLA_Axpyt_h_blk_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &A1, /**/ &A2,
                           b, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &B1, 
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    FLA_Axpyt_internal( FLA_CONJ_TRANSPOSE, alpha, A1, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ A1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* ** */
                                                  B1, 
                              &BB,                B2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

