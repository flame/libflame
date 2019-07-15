/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
  return FLA_QR_UT_internal( A, T,
                             fla_qrut_cntl_leaf );
}

