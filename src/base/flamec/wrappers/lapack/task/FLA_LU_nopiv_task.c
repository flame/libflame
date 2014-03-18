
#include "FLAME.h"

extern fla_lu_t* fla_lu_nopiv_cntl_leaf;

FLA_Error FLA_LU_nopiv_task( FLA_Obj A, fla_lu_t* cntl )
{
  return FLA_LU_nopiv_internal( A,
                                fla_lu_nopiv_cntl_leaf );
}

