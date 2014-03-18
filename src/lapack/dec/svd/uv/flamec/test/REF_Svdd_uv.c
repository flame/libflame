
#include "FLAME.h"

FLA_Error REF_Svdd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  return FLA_Svdd_external( FLA_SVD_VECTORS_ALL, A, s, U, V );
}

