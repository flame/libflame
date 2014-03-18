
#include "FLAME.h"

FLA_Error REF_Hevdr_lv( FLA_Obj A, FLA_Obj l, FLA_Obj Z )
{
  return FLA_Hevdr_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l, Z );
}

