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

FLA_Error FLA_Copyconj_tri_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, void* B_mat )
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

  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
    // lower triangular
    dim_t n_elem_max = m_A;
    switch ( datatype ){

      case FLA_COMPLEX:
      {

        float* buff_B_hip = ( float* ) B_mat;
        float buff_alpha  = -1.0f;

        for ( int j = 0; j < min( n_A, m_A ); j++)
        {
          dim_t n_elem = n_elem_max - j;
          dim_t offset = 2 * ( j * ldim_A + j )  + 1; // plus one for the imaginary part
          rocblas_sscal( handle,         
                         n_elem,
                         &buff_alpha,
                         buff_B_hip + offset, 2 );
        }

        break;
      }
      case FLA_DOUBLE_COMPLEX:
      {

        double* buff_B_hip = ( double* ) B_mat;
        double buff_alpha = -1.0;

        for ( int j = 0; j < min( n_A, m_A); j++)
        {
          dim_t n_elem = n_elem_max - j;
          dim_t offset = 2 * ( j * ldim_A + j ) + 1; // plus one for the imaginary part
          rocblas_dscal( handle,     
                         n_elem,
                         &buff_alpha,
                         buff_B_hip + offset, 2 );
        }

        break;
      }
    }
  }
  else
  {
    // upper triangular
    dim_t n_elem_max = min( m_A, n_A );
    switch ( datatype ){

      case FLA_COMPLEX:
      {

        float* buff_B_hip = ( float* ) B_mat;
        float buff_alpha  = -1.0f;

	for ( int j = 0; j < n_A; j++)
        {
          dim_t n_elem = min( j + 1, n_elem_max );
	  dim_t offset = 2 * ( j * ldim_A ) + 1; // plus one for the imaginary part
	  rocblas_sscal( handle,
                         n_elem,
                         &buff_alpha,
                         buff_B_hip + offset, 2 );
        }

        break;
      }
      case FLA_DOUBLE_COMPLEX:
      {

        double* buff_B_hip = ( double* ) B_mat;
        double buff_alpha = -1.0;

        for ( int j = 0; j < n_A; j++)
        {
          dim_t n_elem = min( j + 1, n_elem_max );
          dim_t offset = 2 * ( j * ldim_A ) + 1; // plus one for the imaginary part
          rocblas_dscal( handle,         
                         n_elem,
                         &buff_alpha,
                         buff_B_hip + offset, 2 );
        }

        break;
      }
    }
  }

  return FLA_SUCCESS;
}

#endif
