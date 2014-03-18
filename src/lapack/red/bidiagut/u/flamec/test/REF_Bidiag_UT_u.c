
#include "FLAME.h"

FLA_Error REF_Bidiag_UT_u( FLA_Obj A, FLA_Obj tu, FLA_Obj tv )
{
  return FLA_Bidiag_blk_external( A, tu, tv );
}

