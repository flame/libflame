
#include "FLAME.h"

FLA_Error REF_Hess( FLA_Obj A, FLA_Obj t, int ilo, int ihi )
{
  return FLA_Hess_blk_external( A, t, ilo, ihi );
}
