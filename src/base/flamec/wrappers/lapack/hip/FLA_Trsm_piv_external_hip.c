/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"

FLA_Error FLA_Trsm_piv_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj p )
{
  FLA_Apply_pivots_unb_external_hip( handle, FLA_LEFT, FLA_NO_TRANSPOSE, 
                                     p, B , B_hip);

  FLA_Trsm_external_hip( handle, FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                         FLA_ONE, A, A_hip, B, B_hip);
  
  return FLA_SUCCESS;
}

#endif
