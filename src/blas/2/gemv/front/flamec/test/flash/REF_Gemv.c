
#include "FLAME.h"

FLA_Error REF_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
  FLA_Gemv_external( transa, alpha, A, x, beta, y );

  return 0;
}

