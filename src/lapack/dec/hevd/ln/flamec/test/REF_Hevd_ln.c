
#include "FLAME.h"

FLA_Error REF_Hevd_ln( FLA_Obj A, FLA_Obj l )
{
  return FLA_Hevd_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, A, l );
}

