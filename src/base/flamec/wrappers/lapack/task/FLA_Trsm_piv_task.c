
#include "FLAME.h"

FLA_Error FLA_Trsm_piv_task( FLA_Obj A, FLA_Obj B, FLA_Obj p, fla_trsm_t* cntl )
{
  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, 
                    p, B );

  FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                     FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                     FLA_ONE, A, B );
  
  return FLA_SUCCESS;
}

