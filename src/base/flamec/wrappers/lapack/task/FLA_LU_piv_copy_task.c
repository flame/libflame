/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LU_piv_copy_task( FLA_Obj A, FLA_Obj p, FLA_Obj U, fla_lu_t* cntl )
{
  FLA_Error r_val;

  r_val = FLA_LU_piv_task( A, p, cntl );

  FLA_Copy_external( A, U );

  return r_val;
}

