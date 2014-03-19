/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_FS_task( FLA_Obj L,
                          FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                                FLA_Obj E, dim_t nb_alg,
                                                           fla_gemm_t* cntl )
{
  FLA_Error info;

  info = FLA_SA_FS_blk( L,
                        D, p, C,
                              E, nb_alg );

  return info;
}

