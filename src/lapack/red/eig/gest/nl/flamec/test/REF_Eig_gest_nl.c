
#include "FLAME.h"

FLA_Error REF_Eig_gest_nl( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
}

