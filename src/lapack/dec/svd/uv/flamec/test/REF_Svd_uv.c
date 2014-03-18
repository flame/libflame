
#include "FLAME.h"

FLA_Error REF_Svd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  return FLA_Svd_external( FLA_SVD_VECTORS_ALL, FLA_SVD_VECTORS_ALL, A, s, U, V );
}

