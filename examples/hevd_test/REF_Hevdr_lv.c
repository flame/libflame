
#include "FLAME.h"

FLA_Error REF_Hevdr_lv( FLA_Obj A, FLA_Obj l, FLA_Obj Z,
                        double* dtime_tred, double* dtime_tevd, double* dtime_appq )

{
  *dtime_tred = 1;
  *dtime_tevd = 1;
  *dtime_appq = 1;

  return FLA_Hevdr_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l, Z );
}

