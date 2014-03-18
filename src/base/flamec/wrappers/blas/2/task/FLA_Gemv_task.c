
#include "FLAME.h"

FLA_Error FLA_Gemv_task( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  return FLA_Gemv_external( transa, alpha, A, x, beta, y );
}

FLA_Error FLA_Gemv_h_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  return FLA_Gemv_external( FLA_CONJ_TRANSPOSE, alpha, A, x, beta, y );
}

FLA_Error FLA_Gemv_n_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  return FLA_Gemv_external( FLA_NO_TRANSPOSE, alpha, A, x, beta, y );
}

FLA_Error FLA_Gemv_t_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  return FLA_Gemv_external( FLA_TRANSPOSE, alpha, A, x, beta, y );
}

