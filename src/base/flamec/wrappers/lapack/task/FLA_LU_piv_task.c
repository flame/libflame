
#include "FLAME.h"

extern fla_lu_t* fla_lu_piv_cntl_leaf;

FLA_Error FLA_LU_piv_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  return FLA_LU_piv_internal( A, p,
                              fla_lu_piv_cntl_leaf );
}

