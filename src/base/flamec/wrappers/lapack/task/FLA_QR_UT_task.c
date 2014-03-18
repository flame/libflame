
#include "FLAME.h"

extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
  return FLA_QR_UT_internal( A, T,
                             fla_qrut_cntl_leaf );
}

