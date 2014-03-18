
#include "FLAME.h"

FLA_Error REF_Trmm( FLA_Side side, FLA_Uplo uploA, FLA_Trans transA, FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{

  FLA_Trmm_external( side, uploA, transA, diagA, alpha, A, B );

  return 0;
}
