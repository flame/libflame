/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include <hip/hip_runtime.h>
#include "rocblas/rocblas.h"

FLA_Error FLA_Copyconj_general_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, void* B_mat )
{

  if ( FLA_Obj_has_zero_dim( A ) )
  {
    return FLA_SUCCESS;
  }

  void* A_mat;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    A_mat = A_hip;
  } 

  dim_t m_A      = FLA_Obj_length( A );
  dim_t n_A      = FLA_Obj_width( A );
  dim_t ldim_A   = FLA_Obj_col_stride( A );

  dim_t elem_size = FLA_Obj_elem_size( A );

  hipStream_t stream;
  rocblas_get_stream( handle, &stream );
  size_t count = elem_size * ldim_A * n_A;
  hipError_t err = hipMemcpyAsync( B_mat,
                                   A_mat,
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

  // early return if for whatever reason a real object was passed
  if ( FLA_Obj_is_real( A ) ) return FLA_SUCCESS;

  FLA_Datatype datatype = FLA_Obj_datatype( A );

  if ( m_A == ldim_A )
  {
    // we can do a fast path that only does one simple rocblas_[s,d]scal invocation
    switch ( datatype ){

      case FLA_COMPLEX:
      {

        float* buff_B_hip = ( float* ) B_mat;
        float buff_alpha  = -1.0f; 

        rocblas_sscal( handle,
                       m_A * n_A,
                       &buff_alpha,
                       buff_B_hip + 1, 2 );

        break;
      }
      case FLA_DOUBLE_COMPLEX:
      {

        double* buff_B_hip = ( double* ) B_mat;
	double buff_alpha = -1.0;

        rocblas_dscal( handle,
                       m_A * n_A,
                       &buff_alpha, 
                       buff_B_hip + 1, 2 );

        break;
      }
    }
  }
  else
  {
    switch ( datatype ){

      case FLA_COMPLEX:
      {

        float* buff_B_hip = ( float* ) B_mat;
        float buff_alpha  = -1.0f;  
        
        rocblas_sscal_strided_batched( handle,
                                       m_A,
                                       &buff_alpha, 
                                       buff_B_hip + 1,
				       2,
				       ldim_A,
				       n_A );

        break;
      }
      case FLA_DOUBLE_COMPLEX:
      {

        double* buff_B_hip = ( double* ) B_mat;
        double buff_alpha = -1.0;

        rocblas_dscal_strided_batched( handle,
                                       m_A,
                                       &buff_alpha,
                                       buff_B_hip + 1,
                                       2,
                                       ldim_A,
                                       n_A );
        break;
      }
    }
  }

  return FLA_SUCCESS;
}

#endif
