/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Hess( FLA_Obj A, FLA_Obj t, integer ilo, integer ihi )
{
  return FLA_Hess_blk_external( A, t, ilo, ihi );
}
