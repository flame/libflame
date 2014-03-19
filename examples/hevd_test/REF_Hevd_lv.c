/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Hevd_lv( FLA_Obj A, FLA_Obj l,
                       double* dtime_tred, double* dtime_tevd, double* dtime_appq )

{
  *dtime_tred = 1;
  *dtime_tevd = 1;
  *dtime_appq = 1;

  return FLA_Hevd_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l );
}

