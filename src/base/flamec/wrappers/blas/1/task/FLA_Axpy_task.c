
#include "FLAME.h"

FLA_Error FLA_Axpy_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpy_t* cntl )
{
  return FLA_Axpy_external( alpha, A, B );
}

