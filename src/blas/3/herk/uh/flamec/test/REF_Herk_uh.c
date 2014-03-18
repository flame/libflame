
#include "FLAME.h"

FLA_Error REF_Herk_uh( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{

  FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, alpha, A, beta, C );
  
  return 0;
}

