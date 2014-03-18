
#include "FLAME.h"

extern fla_chol_t* fla_chol_cntl_leaf;

FLA_Error FLA_Chol_task( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl )
{
  return FLA_Chol_internal( uplo, A,
                            fla_chol_cntl_leaf );
}

FLA_Error FLA_Chol_l_task( FLA_Obj A, fla_chol_t* cntl )
{
  //return FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, A );
  return FLA_Chol_internal( FLA_LOWER_TRIANGULAR, A,
                            fla_chol_cntl_leaf );
}

FLA_Error FLA_Chol_u_task( FLA_Obj A, fla_chol_t* cntl )
{
  //return FLA_Chol_unb_external( FLA_UPPER_TRIANGULAR, A );
  return FLA_Chol_internal( FLA_UPPER_TRIANGULAR, A,
                            fla_chol_cntl_leaf );
}

