/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Trsm_task( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( side, uplo, trans, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_llc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_llh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Trsm_lln_task_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external_ts( FLA_cntl_init_i, FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}
#endif

FLA_Error FLA_Trsm_lln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_llt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_luc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_luh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Trsm_lun_task_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external_ts( FLA_cntl_init_i, FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}
#endif
FLA_Error FLA_Trsm_lun_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_lut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_rlc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_rlh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_rln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_rlt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_ruc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_ruh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_run_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trsm_rut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  return FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

