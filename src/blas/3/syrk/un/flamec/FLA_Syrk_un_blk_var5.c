
#include "FLAME.h"

FLA_Error FLA_Syrk_un_blk_var5( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  dim_t b;

  FLA_Scalr_internal( FLA_UPPER_TRIANGULAR, beta, C,
                      FLA_Cntl_sub_scalr( cntl ) );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* C = C + A1 * A1' */
    FLA_Syrk_internal( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                       alpha, A1, FLA_ONE, C,
                       FLA_Cntl_sub_syrk( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
