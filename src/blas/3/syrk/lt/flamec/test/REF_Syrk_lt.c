
#include "FLAME.h"

FLA_Error REF_Syrk_lt( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{

  FLA_Syrk_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, alpha, A, beta, C );
  
  return 0;
}

