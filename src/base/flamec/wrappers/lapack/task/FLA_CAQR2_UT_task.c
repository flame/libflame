
#include "FLAME.h"

extern fla_caqr2ut_t* fla_caqr2ut_cntl_leaf;

FLA_Error FLA_CAQR2_UT_task( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl )
{
  return FLA_CAQR2_UT_internal( B, D, T,
                                fla_caqr2ut_cntl_leaf );
}

