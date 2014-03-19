/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_uddateut_t* fla_uddateut_cntl_leaf;

FLA_Error FLA_UDdate_UT_task( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl )
{
  return FLA_UDdate_UT_internal( R, C, D, T,
                                 fla_uddateut_cntl_leaf );
}

