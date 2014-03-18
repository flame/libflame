
#include "FLAME.h"

FLA_Error FLA_Trsm_rlt_blk_var3_ht( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Trsm_t* cntl )
{
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  int b;

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ){

    b = 1;

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* ** */
                                              &B1, 
                           BB,                &B2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /* B1 = B1 / tril( A'); */
    FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diagA,
                       alpha, *FLASH_OBJ_PTR_AT( A ), *FLASH_OBJ_PTR_AT( B1 ),
                       NULL );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  B1, 
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

