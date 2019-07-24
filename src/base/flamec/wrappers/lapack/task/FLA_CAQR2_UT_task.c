/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_caqr2ut_t* fla_caqr2ut_cntl_leaf;

FLA_Error FLA_CAQR2_UT_task( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl )
{
  return FLA_CAQR2_UT_internal( B, D, T,
                                fla_caqr2ut_cntl_leaf );
}

