
#include "FLAME.h"

FLA_Error REF_Chol_l( FLA_Obj A )
{
  return FLA_Chol_blk_external( FLA_LOWER_TRIANGULAR, A );
}

