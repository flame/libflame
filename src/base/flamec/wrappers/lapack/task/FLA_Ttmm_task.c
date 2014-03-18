
#include "FLAME.h"

extern fla_ttmm_t* fla_ttmm_cntl_leaf;

FLA_Error FLA_Ttmm_task( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl )
{
  return FLA_Ttmm_internal( uplo, A,
                            fla_ttmm_cntl_leaf );
}

FLA_Error FLA_Ttmm_l_task( FLA_Obj A, fla_ttmm_t* cntl )
{
  //return FLA_Ttmm_unb_external( FLA_LOWER_TRIANGULAR, A );
  return FLA_Ttmm_internal( FLA_LOWER_TRIANGULAR, A,
                            fla_ttmm_cntl_leaf );
}

FLA_Error FLA_Ttmm_u_task( FLA_Obj A, fla_ttmm_t* cntl )
{
  //return FLA_Ttmm_unb_external( FLA_UPPER_TRIANGULAR, A );
  return FLA_Ttmm_internal( FLA_UPPER_TRIANGULAR, A,
                            fla_ttmm_cntl_leaf );
}

