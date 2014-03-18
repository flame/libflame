
#include "FLAME.h"

FLA_Error REF_Bsvd_v( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V )
{
  return FLA_Bsvd_external( uplo, d, e, U, V );
}

