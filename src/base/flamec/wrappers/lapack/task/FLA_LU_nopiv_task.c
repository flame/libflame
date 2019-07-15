/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_lu_t* fla_lu_nopiv_cntl_leaf;

FLA_Error FLA_LU_nopiv_task( FLA_Obj A, fla_lu_t* cntl )
{
  return FLA_LU_nopiv_internal( A,
                                fla_lu_nopiv_cntl_leaf );
}

