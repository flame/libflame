
#include "FLAME.h"

FLA_Error REF_Chol_u( FLA_Obj A )
{
  return FLA_Chol_blk_external( FLA_UPPER_TRIANGULAR, A );
}

