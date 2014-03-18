
#include "FLAME.h"

FLA_Error REF_Sylv_nh( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
}

