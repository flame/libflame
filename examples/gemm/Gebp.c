#include "FLAME.h"
#include "Gemm_prototypes.h"

int Gebp( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, A, B, FLA_ONE, C );

  return FLA_SUCCESS;
}

