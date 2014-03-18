
#include "FLAME.h"

FLA_Error FLA_Copy_task( FLA_Obj A, FLA_Obj B, fla_copy_t* cntl )
{
  return FLA_Copy_external( A, B );
}

