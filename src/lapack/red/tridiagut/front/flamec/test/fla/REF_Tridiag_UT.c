
#include "FLAME.h"

FLA_Error REF_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  return FLA_Tridiag_blk_external( uplo, A, t );
}

