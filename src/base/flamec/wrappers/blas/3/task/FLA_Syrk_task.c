/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Syrk_task( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  return FLA_Syrk_external( uplo, trans, alpha, A, beta, C );
}

FLA_Error FLA_Syrk_ln_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  return FLA_Syrk_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, alpha, A, beta, C );
}

FLA_Error FLA_Syrk_lt_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  return FLA_Syrk_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, alpha, A, beta, C );
}

FLA_Error FLA_Syrk_un_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  return FLA_Syrk_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, alpha, A, beta, C );
}

FLA_Error FLA_Syrk_ut_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
  return FLA_Syrk_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, alpha, A, beta, C );
}

