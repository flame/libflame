/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"

FLA_Error FLA_LU_piv_copy_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p, FLA_Obj U, void* U_hip )
{
  FLA_Error r_val;

  r_val = FLA_LU_piv_blk_external_hip( handle, A, A_hip, p );
  if ( r_val != FLA_SUCCESS )
  {
    return r_val;
  }

  r_val = FLA_Copy_external_hip( handle, A, A_hip, U, U_hip );

  return r_val;
}

#endif
