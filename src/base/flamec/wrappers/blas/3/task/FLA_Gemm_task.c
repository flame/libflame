
#include "FLAME.h"

FLA_Error FLA_Gemm_task( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( transa, transb, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_cc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ch_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_cn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ct_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ht_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nt_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_th_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tt_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

