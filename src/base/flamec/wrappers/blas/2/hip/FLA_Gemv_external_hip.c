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

FLA_Error FLA_Gemv_external_hip( rocblas_handle handle, FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj x, void* x_hip, FLA_Obj beta, FLA_Obj y, void* y_hip )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          ldim_A;
  int          inc_x;
  int          inc_y;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemv_check( transa, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_col_stride( A );

  inc_x    = 1;
  inc_y    = 1;

  void* A_mat = NULL;
  void* x_vec = NULL;
  void* y_vec = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    x_vec = FLA_Obj_buffer_at_view( x );
    y_vec = FLA_Obj_buffer_at_view( y );
  }
  else
  {
    A_mat = A_hip;
    x_vec = x_hip;
    y_vec = y_hip;
  }

  FLA_Trans trans_a_corr = transa;
  FLA_Bool conj_no_trans_a = FALSE;

  void* A_mat_corr = NULL;
  if ( FLA_Obj_is_complex( A ) && transa == FLA_CONJ_NO_TRANSPOSE )
  {
    // must correct by copying to temporary buffer and conjugating there
    trans_a_corr = FLA_NO_TRANSPOSE;
    conj_no_trans_a = TRUE;

    dim_t elem_size = FLA_Obj_elem_size( A );
    size_t count = elem_size * ldim_A * n_A;
    hipMalloc( &A_mat_corr, count );
    FLA_Copyconj_general_external_hip( handle, A, A_hip, A_mat_corr );
    A_mat = A_mat_corr;
  }

  rocblas_operation blas_transa = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );

  switch( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    rocblas_sgemv( handle,
		   blas_transa,
                   m_A,
                   n_A, 
                   buff_alpha,                   
                   ( float * ) A_mat, ldim_A,
                   ( float * ) x_vec, inc_x,
                   buff_beta,
                   ( float * ) y_vec, inc_y );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    rocblas_dgemv( handle,
		   blas_transa,
                   m_A,
                   n_A, 
                   buff_alpha,                   
                   ( double * ) A_mat, ldim_A,
                   ( double * ) x_vec, inc_x,
                   buff_beta,
                   ( double * ) y_vec, inc_y );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_alpha = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex *buff_beta  = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( beta );

    rocblas_cgemv( handle,
		   blas_transa,
                   m_A,
                   n_A, 
                   buff_alpha,                   
                   ( rocblas_float_complex * ) A_mat, ldim_A,
                   ( rocblas_float_complex * ) x_vec, inc_x,
                   buff_beta,
                   ( rocblas_float_complex * ) y_vec, inc_y );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_alpha = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex *buff_beta  = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    rocblas_zgemv( handle,
		   blas_transa,
                   m_A,
                   n_A, 
                   buff_alpha,                   
                   ( rocblas_double_complex * ) A_mat, ldim_A,
                   ( rocblas_double_complex * ) x_vec, inc_x,
                   buff_beta,
                   ( rocblas_double_complex * ) y_vec, inc_y );

    break;
  }

  }

  if( conj_no_trans_a )
  {
    hipFree( A_mat_corr );
  }

  return FLA_SUCCESS;
}

#endif
