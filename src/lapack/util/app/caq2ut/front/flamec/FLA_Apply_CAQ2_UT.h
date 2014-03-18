
#include "FLA_Apply_CAQ2_UT_lhfc.h"

FLA_Error FLA_Apply_CAQ2_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                      FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                       FLA_Obj E,
                                      fla_apcaq2ut_t* cntl );

FLA_Error FLA_Apply_CAQ2_UT_lhfc( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                   FLA_Obj E,
                                  fla_apcaq2ut_t* cntl );
