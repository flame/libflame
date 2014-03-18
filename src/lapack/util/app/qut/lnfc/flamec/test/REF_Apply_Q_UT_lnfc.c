
#include "FLAME.h"

FLA_Error REF_Apply_Q_UT_lnfc( FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  return FLA_Apply_Q_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, A, t, B );
}

