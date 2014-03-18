
#include "FLAME.h"

FLA_Error REF_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Axpyt_external( trans, alpha, A, B );

  return 0;
}

