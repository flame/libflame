/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Trsv_task( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( uplo, trans, diag, A, x );
}

FLA_Error FLA_Trsv_lc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, A, x );
}

FLA_Error FLA_Trsv_ln_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, A, x );
}

FLA_Error FLA_Trsv_lt_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, A, x );
}

FLA_Error FLA_Trsv_uc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, A, x );
}

FLA_Error FLA_Trsv_un_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, A, x );
}

FLA_Error FLA_Trsv_ut_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
  return FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, A, x );
}

