
#include "FLAME.h"

FLA_Error REF_Scalr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  return FLA_Scalr_external( uplo, alpha, A );
}

