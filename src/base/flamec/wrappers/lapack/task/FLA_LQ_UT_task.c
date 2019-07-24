/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_lqut_t* fla_lqut_cntl_leaf;

FLA_Error FLA_LQ_UT_task( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl )
{
  return FLA_LQ_UT_internal( A, T,
                             fla_lqut_cntl_leaf );
}

