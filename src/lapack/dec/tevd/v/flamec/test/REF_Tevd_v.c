
#include "FLAME.h"

FLA_Error REF_Tevd_v( FLA_Obj d, FLA_Obj e, FLA_Obj U )
{
  return FLA_Tevd_external( FLA_EVD_WITH_VECTORS, d, e, U );
}

