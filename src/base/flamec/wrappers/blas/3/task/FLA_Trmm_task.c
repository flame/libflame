/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Trmm_task( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( side, uplo, trans, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_luc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_luh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lun_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_ruc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_ruh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_run_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

