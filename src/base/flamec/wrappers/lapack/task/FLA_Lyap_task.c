
#include "FLAME.h"

extern fla_lyap_t* fla_lyap_cntl_leaf;

FLA_Error FLA_Lyap_task( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
  return FLA_Lyap_internal( trans, isgn, A, C, scale,
                            fla_lyap_cntl_leaf );
}

FLA_Error FLA_Lyap_n_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
  return FLA_Lyap_internal( FLA_NO_TRANSPOSE, isgn, A, C, scale,
                            fla_lyap_cntl_leaf );
}

FLA_Error FLA_Lyap_h_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
  return FLA_Lyap_internal( FLA_CONJ_TRANSPOSE, isgn, A, C, scale,
                            fla_lyap_cntl_leaf );
}
