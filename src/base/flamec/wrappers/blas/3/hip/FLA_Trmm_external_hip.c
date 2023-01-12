/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include <hip/hip_runtime.h>
#include "rocblas/rocblas.h"

FLA_Error FLA_Trmm_external_hip( rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{
  FLA_Datatype datatype;
  int          n_A;
  int          m_B, n_B;
  int          ldim_A;
  int          ldim_B;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trmm_check( side, uplo, trans, diag, alpha, A, B );

  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_col_stride( A );
  n_A      = FLA_Obj_width( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_col_stride( B );

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
    hipStream_t stream;
    rocblas_get_stream( handle, &stream );
    hipMallocAsync( &A_mat_corr, count, stream );
    FLA_Copyconj_tri_external_hip( handle, uplo, A, A_hip, A_mat_corr );
    A_mat = A_mat_corr;
  }

  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );
  rocblas_side blas_side = FLA_Param_map_flame_to_rocblas_side( side );
  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );
  rocblas_diagonal blas_diag = FLA_Param_map_flame_to_rocblas_diag( diag );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    rocblas_strmm( handle,
                   blas_side,
                   blas_uplo, 
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( float * ) A_mat, ldim_A,
                   ( float * ) B_mat, ldim_B );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    rocblas_dtrmm( handle,
                   blas_side,
                   blas_uplo, 
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( double * ) A_mat, ldim_A,
                   ( double * ) B_mat, ldim_B );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_alpha = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( alpha );

    rocblas_ctrmm( handle,
                   blas_side,
                   blas_uplo, 
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( rocblas_float_complex * ) A_mat, ldim_A,
                   ( rocblas_float_complex * ) B_mat, ldim_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_alpha = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    rocblas_ztrmm( handle,
                   blas_side,
                   blas_uplo, 
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( rocblas_double_complex * ) A_mat, ldim_A,
                   ( rocblas_double_complex * ) B_mat, ldim_B );

    break;
  }

  }

  if( conj_no_trans_a )
  {
    hipStream_t stream;
    rocblas_get_stream( handle, &stream );
    hipFreeAsync( stream, A_mat_corr );
  }

  return FLA_SUCCESS;
}

#endif
