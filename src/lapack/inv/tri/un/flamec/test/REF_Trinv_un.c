
#include "FLAME.h"

FLA_Error REF_Trinv_un( FLA_Obj A )
{
  return FLA_Trinv_blk_external( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
}

