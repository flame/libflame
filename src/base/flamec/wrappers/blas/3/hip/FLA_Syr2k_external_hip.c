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

FLA_Error FLA_Syr2k_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip )
{
  FLA_Datatype datatype;
  int          k_AB;
  int          m_A, n_A;
  int          m_C;
  int          ldim_A;
  int          ldim_B;
  int          ldim_C;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Syr2k_check( uplo, trans, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_col_stride( A );

  ldim_B   = FLA_Obj_col_stride( B );

  m_C      = FLA_Obj_length( C );
  ldim_C   = FLA_Obj_col_stride( C );

  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
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

  FLA_Trans trans_a_corr = trans;
  FLA_Bool conj_no_trans_a = FALSE;

  void* A_mat_corr = NULL;
  if ( FLA_Obj_is_complex( A ) && trans == FLA_CONJ_NO_TRANSPOSE )
  {
    // must correct by copying to temporary buffer and conjugating there
    trans_a_corr = FLA_NO_TRANSPOSE;
    conj_no_trans_a = TRUE;

    dim_t elem_size = FLA_Obj_elem_size( A );
    size_t count = elem_size * ldim_A * n_A;
    hipMalloc( &A_mat_corr, count );
    FLA_Copyconj_tri_external_hip( handle, uplo, A, A_hip, A_mat_corr );
    A_mat = A_mat_corr;
  }

  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );
  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    rocblas_ssyr2k( handle,
                    blas_uplo,
                    blas_trans,
                    m_C,
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

    rocblas_dsyr2k( handle,
                    blas_uplo,
                    blas_trans,
                    m_C,
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

    rocblas_csyr2k( handle,
                    blas_uplo,
                    blas_trans,
                    m_C,
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

    rocblas_zsyr2k( handle,
                    blas_uplo,
                    blas_trans,
                    m_C,
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

  return FLA_SUCCESS;
}

#endif
