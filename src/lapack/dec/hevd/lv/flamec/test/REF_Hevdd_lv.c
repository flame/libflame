
#include "FLAME.h"

FLA_Error REF_Hevdd_lv( FLA_Obj A, FLA_Obj l )
{
  return FLA_Hevdd_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l );
}

