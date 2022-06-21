/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "rocblas.h"

FLA_Error FLA_Trsv_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_hip, FLA_Obj x, void* x_hip ) 
{
  FLA_Datatype datatype;
  int          m_A;
  int          ldim_A;
  int          inc_x;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  ldim_A   = FLA_Obj_col_stride( A );

  inc_x    = 1;

  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );
  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans, FLA_Obj_is_real( A ) );
  rocblas_diagonal blas_diag = FLA_Param_map_flame_to_rocblas_diag( diag );

  void* A_mat = NULL;
  void* x_vec = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip( ) )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    x_vec = FLA_Obj_buffer_at_view( x );
  }
  else
  {
    A_mat = A_hip;
    x_vec = x_hip;
  }

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

  return FLA_SUCCESS;
}

#endif
