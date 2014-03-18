
#include "FLAME.h"

FLA_Error FLA_Scalr_task( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl )
{
  return FLA_Scalr_external( uplo, alpha, A );
}

FLA_Error FLA_Scalr_l_task( FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl )
{
  return FLA_Scalr_external( FLA_LOWER_TRIANGULAR, alpha, A );
}

FLA_Error FLA_Scalr_u_task( FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl )
{
  return FLA_Scalr_external( FLA_UPPER_TRIANGULAR, alpha, A );
}

