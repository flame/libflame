
#include "FLAME.h"

FLA_Error REF_Scal( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Scal_external( alpha, A, B );

  return 0;
}

