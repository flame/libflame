
#include "FLAME.h"

extern fla_qr2ut_t* fla_qr2ut_cntl_leaf;

FLA_Error FLA_QR2_UT_task( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl )
{
  return FLA_QR2_UT_internal( B, D, T,
                              fla_qr2ut_cntl_leaf );
}

