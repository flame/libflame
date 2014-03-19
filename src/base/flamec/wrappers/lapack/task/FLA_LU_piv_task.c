/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_lu_t* fla_lu_piv_cntl_leaf;

FLA_Error FLA_LU_piv_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  return FLA_LU_piv_internal( A, p,
                              fla_lu_piv_cntl_leaf );
}

