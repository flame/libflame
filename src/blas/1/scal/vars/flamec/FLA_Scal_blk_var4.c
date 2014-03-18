
#include "FLAME.h"

FLA_Error FLA_Scal_blk_var4( FLA_Obj alpha, FLA_Obj A, fla_scal_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  dim_t b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &A1, /**/ &A2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    FLA_Scal_internal( alpha, A1,
                       FLA_Cntl_sub_scal( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ A1, A2,
                              FLA_RIGHT );
  }

  return FLA_SUCCESS;
}

