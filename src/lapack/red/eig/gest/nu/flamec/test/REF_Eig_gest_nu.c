
#include "FLAME.h"

FLA_Error REF_Eig_gest_nu( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
}

