
#include "FLAME.h"

FLA_Error FLA_Scal_task( FLA_Obj alpha, FLA_Obj A, fla_scal_t* cntl )
{
  return FLA_Scal_external( alpha, A );
}

