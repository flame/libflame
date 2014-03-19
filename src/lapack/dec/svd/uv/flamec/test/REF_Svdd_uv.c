/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Svdd_uv( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  return FLA_Svdd_external( FLA_SVD_VECTORS_ALL, A, s, U, V );
}

