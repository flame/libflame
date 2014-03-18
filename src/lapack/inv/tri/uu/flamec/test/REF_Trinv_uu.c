
#include "FLAME.h"

FLA_Error REF_Trinv_uu( FLA_Obj A )
{
  return FLA_Trinv_blk_external( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A );
}

