/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_LU_task( FLA_Obj U,
                          FLA_Obj D, FLA_Obj p, FLA_Obj L, dim_t nb_alg,
                                                           fla_lu_t* cntl )
{
  FLA_Error info;

  info = FLA_SA_LU_blk( U,
                        D, p, L, nb_alg );

  return info;
}

