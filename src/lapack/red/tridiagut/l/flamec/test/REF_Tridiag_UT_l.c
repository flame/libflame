/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Tridiag_UT_l( FLA_Obj A, FLA_Obj t )
{
  return FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, A, t );
}

