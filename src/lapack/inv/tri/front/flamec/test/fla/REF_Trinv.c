
#include "FLAME.h"

FLA_Error REF_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
  return FLA_Trinv_blk_external( uplo, diag, A );
}

