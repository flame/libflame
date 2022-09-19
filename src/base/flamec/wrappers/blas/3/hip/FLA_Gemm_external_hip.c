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

FLA_Error FLA_Gemm_external_hip( rocblas_handle handle, FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip )
{
  FLA_Datatype datatype;
  int          k_AB;
  int          m_A, n_A;
  int          n_B;
  int          m_C, n_C;
  int          ldim_A;
  int          ldim_B;
  int          ldim_C;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemm_check( transa, transb, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  if ( FLA_Obj_has_zero_dim( A ) || FLA_Obj_has_zero_dim( B ) )
  {
    FLA_Scal_external_hip( handle, beta, C, C_hip );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_col_stride( A );

  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_col_stride( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  ldim_C   = FLA_Obj_col_stride( C );

  if ( transa == FLA_NO_TRANSPOSE || transa == FLA_CONJ_NO_TRANSPOSE )
    k_AB = n_A;
  else
    k_AB = m_A;

  void* A_mat = NULL;
  void* B_mat = NULL;
  void* C_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    B_mat = FLA_Obj_buffer_at_view( B );
    C_mat = FLA_Obj_buffer_at_view( C );
  }
  else
  {
    A_mat = A_hip;
    B_mat = B_hip;
    C_mat = C_hip;
  }

  FLA_Trans trans_a_corr = transa;
  FLA_Trans trans_b_corr = transb;
  FLA_Bool conj_no_trans_a = FALSE;
  FLA_Bool conj_no_trans_b = FALSE;

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

  void* B_mat_corr = NULL;
  if ( FLA_Obj_is_complex( B ) && transb == FLA_CONJ_NO_TRANSPOSE )
  {
    // must correct by copying to temporary buffer and conjugating there
    trans_b_corr = FLA_NO_TRANSPOSE;
    conj_no_trans_b = TRUE;

    dim_t elem_size = FLA_Obj_elem_size( B );
    size_t count = elem_size * ldim_B * n_B;
    hipMalloc( &B_mat_corr, count );
    FLA_Copyconj_general_external_hip( handle, B, B_hip, B_mat_corr );
    B_mat = B_mat_corr;
  }

  rocblas_operation blas_transa = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );
  rocblas_operation blas_transb = FLA_Param_map_flame_to_rocblas_trans( trans_b_corr, FLA_Obj_is_real( B ) );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    rocblas_sgemm( handle,
                   blas_transa,
                   blas_transb,
                   m_C,
                   n_C,
                   k_AB,
                   buff_alpha,
                   ( float * ) A_mat, ldim_A,
                   ( float * ) B_mat, ldim_B,
                   buff_beta,
                   ( float * ) C_mat, ldim_C );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    rocblas_dgemm( handle,
                   blas_transa,
                   blas_transb,
                   m_C,
                   n_C,
                   k_AB,
                   buff_alpha,
                   ( double * ) A_mat, ldim_A,
                   ( double * ) B_mat, ldim_B,
                   buff_beta,
                   ( double * ) C_mat, ldim_C );
    
    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_alpha = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex *buff_beta  = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( beta );

    rocblas_cgemm( handle,
                   blas_transa,
                   blas_transb,
                   m_C,
                   n_C,
                   k_AB,
                   buff_alpha,
                   ( rocblas_float_complex * ) A_mat, ldim_A,
                   ( rocblas_float_complex * ) B_mat, ldim_B,
                   buff_beta,
                   ( rocblas_float_complex * ) C_mat, ldim_C );
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_alpha = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex *buff_beta  = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    rocblas_zgemm( handle,
                   blas_transa,
                   blas_transb,
                   m_C,
                   n_C,
                   k_AB,
                   buff_alpha,
                   ( rocblas_double_complex * ) A_mat, ldim_A,
                   ( rocblas_double_complex * ) B_mat, ldim_B,
                   buff_beta,
                   ( rocblas_double_complex * ) C_mat, ldim_C );
    
    break;
  }

  }

  if( conj_no_trans_a )
  {
    hipFree( A_mat_corr );
  }
  if( conj_no_trans_b )
  {
    hipFree( B_mat_corr );
  }

  return FLA_SUCCESS;
}

#endif
