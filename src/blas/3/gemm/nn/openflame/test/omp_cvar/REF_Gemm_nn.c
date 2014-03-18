
#include "FLAME.h"

FLA_Error REF_Gemm( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Gemm( transa, transb, alpha, A, B, beta, C );
  
  return 0;
}

