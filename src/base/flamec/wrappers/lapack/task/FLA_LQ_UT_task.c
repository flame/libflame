
#include "FLAME.h"

extern fla_lqut_t* fla_lqut_cntl_leaf;

FLA_Error FLA_LQ_UT_task( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl )
{
  return FLA_LQ_UT_internal( A, T,
                             fla_lqut_cntl_leaf );
}

