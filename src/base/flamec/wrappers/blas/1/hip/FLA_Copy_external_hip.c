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

FLA_Error FLA_Copy_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Copy_check( A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  dim_t elem_size = FLA_Obj_elem_size( B );
  dim_t ldim_A    = FLA_Obj_col_stride( A );

  dim_t m_B       = FLA_Obj_length( B );
  dim_t n_B       = FLA_Obj_width( B );
  dim_t ldim_B    = FLA_Obj_col_stride( B );

  void* A_mat = NULL;
  void* B_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    B_mat = FLA_Obj_buffer_at_view( B );
  }
  else
  {
    A_mat = A_hip;
    B_mat = B_hip;
  }

  hipStream_t stream;
  rocblas_get_stream( handle, &stream );

  if ( ldim_A == ldim_B && m_B == ldim_A )
  {
    // clean setup - just transfer the entire block
    size_t count = elem_size * ldim_A * n_B;
    hipError_t err = hipMemcpyAsync( B_mat,
                                     A_mat,
                                     count,
                                     hipMemcpyDeviceToDevice,
                                     stream );
    if ( err != hipSuccess )
    {
      fprintf( stderr,
               "Failure to memcpy D2D in HIP. Size=%ld, err=%d\n",
               count, err );
      return FLA_FAILURE;
    }

  }
  else
  {

    size_t dpitch   = elem_size * ldim_B;
    size_t spitch   = elem_size * ldim_A;
    size_t width    = elem_size * m_B;
    size_t height   = n_B;

    // hipMemcpy2D assumes row-major layout
    hipError_t err = hipMemcpy2DAsync( B_mat,
                                       dpitch,
                                       A_mat,
                                       spitch,
                                       width,
                                       height,
                                       hipMemcpyDeviceToDevice,
                                       stream );

    if ( err != hipSuccess )
    {
      fprintf( stderr, "Failure to 2D memcpy D2D in HIP: %d\n", err);
      return FLA_FAILURE;
    }
  }

  return FLA_SUCCESS;
}

#endif
