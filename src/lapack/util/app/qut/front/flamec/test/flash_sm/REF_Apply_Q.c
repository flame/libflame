
#include "FLAME.h"

FLA_Error REF_Apply_Q( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  return FLA_Apply_Q_blk_external( side, uplo, trans, A, t, B );
}

