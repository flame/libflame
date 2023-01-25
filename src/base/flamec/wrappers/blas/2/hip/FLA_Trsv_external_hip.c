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

FLA_Error FLA_Trsv_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_hip, FLA_Obj x, void* x_hip ) 
{
  FLA_Datatype datatype;
  int          n_A;
  int          m_A;
  int          ldim_A;
  int          inc_x;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  m_A      = FLA_Obj_length( A );
  ldim_A   = FLA_Obj_col_stride( A );

  inc_x    = 1;

  void* A_mat = NULL;
  void* x_vec = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    x_vec = FLA_Obj_buffer_at_view( x );
  }
  else
  {
    A_mat = A_hip;
    x_vec = x_hip;
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

  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );
  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans_a_corr, FLA_Obj_is_real( A ) );
  rocblas_diagonal blas_diag = FLA_Param_map_flame_to_rocblas_diag( diag );

  switch( datatype ){

  case FLA_FLOAT:
  {
    rocblas_strsv( handle, blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_A,
                   ( float * ) A_mat, ldim_A,
                   ( float * ) x_vec, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    rocblas_dtrsv( handle, blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_A,
                   ( double * ) A_mat, ldim_A,
                   ( double * ) x_vec, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_ctrsv( handle, blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_A,
                   ( rocblas_float_complex * ) A_mat, ldim_A,
                   ( rocblas_float_complex * ) x_vec, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_ztrsv( handle, blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_A,
                   ( rocblas_double_complex * ) A_mat, ldim_A,
                   ( rocblas_double_complex * ) x_vec, inc_x );

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
