/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"
#include "Gemm_prototypes.h"

int Gebp( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, A, B, FLA_ONE, C );

  return FLA_SUCCESS;
}

