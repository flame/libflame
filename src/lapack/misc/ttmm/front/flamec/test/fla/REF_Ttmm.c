
#include "FLAME.h"

FLA_Error REF_Ttmm( FLA_Uplo uplo, FLA_Obj A )
{
  return FLA_Ttmm_blk_external( uplo, A );
}

