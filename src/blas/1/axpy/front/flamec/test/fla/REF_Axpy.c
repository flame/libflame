
#include "FLAME.h"

FLA_Error REF_Axpy( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  return FLA_Axpy_external( alpha, A, B );
}

