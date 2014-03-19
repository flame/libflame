/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Apply_Q_UT_lnfc( FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  return FLA_Apply_Q_blk_external( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, A, t, B );
}

