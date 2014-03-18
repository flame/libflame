
#include "FLAME.h"

FLA_Error REF_Tridiag_UT_l( FLA_Obj A, FLA_Obj t )
{
  return FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, A, t );
}

