/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Symm_task( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
  return FLA_Symm_external( side, uplo, alpha, A, B, beta, C );
}

FLA_Error FLA_Symm_ll_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
  return FLA_Symm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, alpha, A, B, beta, C );
}

FLA_Error FLA_Symm_lu_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
  return FLA_Symm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, alpha, A, B, beta, C );
}

FLA_Error FLA_Symm_rl_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
  return FLA_Symm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, alpha, A, B, beta, C );
}

FLA_Error FLA_Symm_ru_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
  return FLA_Symm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, alpha, A, B, beta, C );
}

