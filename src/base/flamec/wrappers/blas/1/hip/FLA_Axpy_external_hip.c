/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "rocblas/rocblas.h"

FLA_Error FLA_Axpy_external_hip( rocblas_handle handle, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Axpy_check( alpha, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  FLA_Datatype datatype = FLA_Obj_datatype( A );

  rocblas_stride ldim_A = FLA_Obj_col_stride( A );
  rocblas_int inc_A     = 1;

  rocblas_int m_B       = FLA_Obj_length( B );
  rocblas_int n_B       = FLA_Obj_width( B );
  rocblas_stride ldim_B = FLA_Obj_col_stride( B );
  rocblas_int inc_B     = 1;

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

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_alpha = ( float* ) FLA_FLOAT_PTR( alpha );
    float* buff_A_hip = ( float* ) A_mat;
    float* buff_B_hip = ( float* ) B_mat;

    rocblas_status err = rocblas_saxpy_strided_batched( handle,
                                   m_B,
                                   buff_alpha,
                                   buff_A_hip,
                                   inc_A,
                                   ldim_A,
                                   buff_B_hip,
                                   inc_B,
                                   ldim_B,
                                   n_B );

    if ( err != rocblas_status_success )
    {
      fprintf( stderr, "Failure in rocblas_saxpy_strided_batched: %d",
               err );
      return FLA_FAILURE;
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_alpha = ( double* ) FLA_DOUBLE_PTR( alpha );
    double* buff_A_hip = ( double* ) A_mat;
    double* buff_B_hip = ( double* ) B_mat;

    rocblas_status err = rocblas_daxpy_strided_batched( handle,
                                   m_B,
                                   buff_alpha,
                                   buff_A_hip,
                                   inc_A,
                                   ldim_A,
                                   buff_B_hip,
                                   inc_B,
                                   ldim_B,
                                   n_B );

    if ( err != rocblas_status_success )
    {
      fprintf( stderr, "Failure in rocblas_daxpy_strided_batched: %d",
               err );
      return FLA_FAILURE;
    }

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_alpha = ( rocblas_float_complex* ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex* buff_A_hip = ( rocblas_float_complex* ) A_mat;
    rocblas_float_complex* buff_B_hip = ( rocblas_float_complex* ) B_mat;

    rocblas_status err = rocblas_caxpy_strided_batched( handle,
                                   m_B,
                                   buff_alpha,
                                   buff_A_hip,
                                   inc_A,
                                   ldim_A,
                                   buff_B_hip,
                                   inc_B,
                                   ldim_B,
                                   n_B );

    if ( err != rocblas_status_success )
    {
      fprintf( stderr, "Failure in rocblas_caxpy_strided_batched: %d",
               err );
      return FLA_FAILURE;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_alpha = ( rocblas_double_complex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex* buff_A_hip = ( rocblas_double_complex* ) A_mat;
    rocblas_double_complex* buff_B_hip = ( rocblas_double_complex* ) B_mat;

    rocblas_status err = rocblas_zaxpy_strided_batched( handle,
                                   m_B,
                                   buff_alpha,
                                   buff_A_hip,
                                   inc_A,
                                   ldim_A,
                                   buff_B_hip,
                                   inc_B,
                                   ldim_B,
                                   n_B );

    if ( err != rocblas_status_success )
    {
      fprintf( stderr, "Failure in rocblas_zaxpy_strided_batched: %d",
               err );
      return FLA_FAILURE;
    }

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
