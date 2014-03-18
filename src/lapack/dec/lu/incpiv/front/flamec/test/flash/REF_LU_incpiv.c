
#include "FLAME.h"

FLA_Error REF_LU_incpiv( FLA_Obj A, FLA_Obj p )
{
  return FLA_LU_piv_blk_external( A, p );
}

