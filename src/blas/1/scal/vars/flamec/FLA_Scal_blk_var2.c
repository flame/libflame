
#include "FLAME.h"

FLA_Error FLA_Scal_blk_var2( FLA_Obj alpha, FLA_Obj A, fla_scal_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    FLA_Scal_internal( alpha, A1,
                       FLA_Cntl_sub_scal( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

