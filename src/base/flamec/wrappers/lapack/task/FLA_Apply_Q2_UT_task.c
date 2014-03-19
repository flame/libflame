/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_apq2ut_t* fla_apq2ut_cntl_leaf;

FLA_Error FLA_Apply_Q2_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                 FLA_Obj E, fla_apq2ut_t* cntl )
{
  return FLA_Apply_Q2_UT_internal( side, trans, direct, storev,
                                   D, T, W, C, E,
                                   fla_apq2ut_cntl_leaf );
}

FLA_Error FLA_Apply_Q2_UT_lhfc_task( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                      FLA_Obj E, fla_apq2ut_t* cntl )
{
  return FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                   D, T, W, C, E,
                                   fla_apq2ut_cntl_leaf );
}

