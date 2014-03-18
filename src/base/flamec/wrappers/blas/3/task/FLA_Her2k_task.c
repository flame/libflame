
#include "FLAME.h"

FLA_Error FLA_Her2k_task( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
  return FLA_Her2k_external( uplo, trans, alpha, A, B, beta, C );
}

FLA_Error FLA_Her2k_ln_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
  return FLA_Her2k_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Her2k_lh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
  return FLA_Her2k_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Her2k_un_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
  return FLA_Her2k_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Her2k_uh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
  return FLA_Her2k_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

