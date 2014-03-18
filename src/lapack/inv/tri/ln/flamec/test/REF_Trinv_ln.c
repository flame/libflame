
#include "FLAME.h"

FLA_Error REF_Trinv_ln( FLA_Obj A )
{
  //return FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
  return FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
}

