
#include "FLAME.h"

FLA_Error REF_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_blk_external( transa, transb, isgn, A, B, C, scale );
}

