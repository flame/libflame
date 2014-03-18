
#include "FLAME.h"

FLA_Error REF_Trinv_lu( FLA_Obj A )
{
  return FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_UNIT_DIAG, A );
}

