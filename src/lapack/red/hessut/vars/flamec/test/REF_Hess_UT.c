
#include "FLAME.h"

FLA_Error REF_Hess_UT( FLA_Obj A, FLA_Obj t )
{
  int ilo = 0;
  int ihi = FLA_Obj_length( A ) - 1;

  return FLA_Hess_blk_external( A, t, ilo, ihi );
}

