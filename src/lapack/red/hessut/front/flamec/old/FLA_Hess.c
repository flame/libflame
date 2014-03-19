/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_hess_t*      fla_hess_cntl;
extern fla_blocksize_t* fla_hess_blocksize;

FLA_Error FLA_Hess( FLA_Obj A, FLA_Obj t, int ilo, int ihi )
{
  FLA_Error r_val;

  // Note the following constraints:  0 <= ilo <= n-1,  0 <= ihi <= n-1.
  r_val = FLA_Hess_internal( A, t, ilo, ihi, fla_hess_cntl );

  return r_val;
}

