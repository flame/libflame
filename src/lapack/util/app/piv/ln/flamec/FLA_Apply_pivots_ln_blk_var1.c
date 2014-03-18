
#include "FLAME.h"

FLA_Error FLA_Apply_pivots_ln_blk_var1( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
  FLA_Obj AL,  AR,       A0,  A1,  A2;

  dim_t b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ) {

    b = FLA_Determine_blocksize( AR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* Apply pivots to each column panel */
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p, A1,
                               FLA_Cntl_sub_appiv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}
