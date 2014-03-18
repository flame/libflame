
#include "FLAME.h"

extern fla_tpose_t* fla_tpose_cntl;

FLA_Error FLA_Transpose( FLA_Obj A )
{
  FLA_Error r_val;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Transpose_check( A );

  r_val = FLA_Transpose_blk_var2( A, fla_tpose_cntl );

  return r_val;
}

