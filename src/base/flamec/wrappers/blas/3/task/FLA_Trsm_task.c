
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

