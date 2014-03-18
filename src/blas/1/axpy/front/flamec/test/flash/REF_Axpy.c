
#include "FLAME.h"

FLA_Error REF_Axpy( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Axpy_external( alpha, A, B );

  return 0;
}

