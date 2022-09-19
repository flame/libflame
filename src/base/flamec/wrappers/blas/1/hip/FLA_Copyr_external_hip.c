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

FLA_Error FLA_Copyr_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{
  dim_t          m_A, n_A;
  dim_t          cs_A;
  dim_t          cs_B;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Copyr_check( uplo, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  cs_B     = FLA_Obj_col_stride( B );

  dim_t elem_size = FLA_Obj_elem_size( B );

  void* A_mat;
  void* B_mat;
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

  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
    // lower triangular
    dim_t n_elem_max = m_A;
    for ( int j = 0; j < min( n_A, m_A ); j++)
    {
      dim_t n_elem = n_elem_max - j;
      dim_t count = elem_size * n_elem;
      dim_t offset_A = elem_size * ( j * cs_A + j );
      dim_t offset_B = elem_size * ( j * cs_B + j );
      hipError_t err = hipMemcpyAsync( B_mat + offset_B,
                                       A_mat + offset_A,
                                       count,
                                       hipMemcpyDeviceToDevice,
                                       stream );
      if ( err != hipSuccess )
      {
        fprintf( stderr,
                 "Failure to copy block HIP2HIP device. Size=%ld, err=%d\n",
                 count, err );
        return FLA_FAILURE;
      }
    }
  }
  else
  {
    // upper triangular
    dim_t n_elem_max = min( m_A, n_A );
    for ( int j = 0; j < n_A; j++)
    {
      dim_t n_elem = min( j + 1, n_elem_max );
      dim_t count = elem_size * n_elem;
      dim_t offset_A = elem_size * ( j * cs_A );
      dim_t offset_B = elem_size * ( j * cs_B );
      hipError_t err = hipMemcpyAsync( B_mat + offset_B,
                                       A_mat + offset_A,
                                       count,
                                       hipMemcpyDeviceToDevice,
                                       stream );
      if ( err != hipSuccess )
      {
        fprintf( stderr,
                 "Failure to copy block HIP2HIP device. Size=%ld, err=%d\n",
                 count, err );
        return FLA_FAILURE;
      }
    }
  }
  
  return FLA_SUCCESS;
}

#endif
