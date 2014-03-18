
#include "FLAME.h"

extern fla_apcaq2ut_t* fla_apcaq2ut_cntl_leaf;

FLA_Error FLA_Apply_CAQ2_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                   FLA_Obj E, fla_apcaq2ut_t* cntl )
{
  return FLA_Apply_CAQ2_UT_internal( side, trans, direct, storev,
                                     D, T, W, C, E,
                                     fla_apcaq2ut_cntl_leaf );
}

FLA_Error FLA_Apply_CAQ2_UT_lhfc_task( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                        FLA_Obj E, fla_apcaq2ut_t* cntl )
{
  return FLA_Apply_CAQ2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                     D, T, W, C, E,
                                     fla_apcaq2ut_cntl_leaf );
}

