
#include "FLAME.h"

FLA_Error FLA_Copyr_task( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( uplo, A, B );
}

FLA_Error FLA_Copyr_l_task( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Copyr_u_task( FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
  return FLA_Copyr_external( FLA_UPPER_TRIANGULAR, A, B );
}

