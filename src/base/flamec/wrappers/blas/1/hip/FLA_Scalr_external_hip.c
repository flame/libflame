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

FLA_Error FLA_Scalr_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_hip )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          ldim_A, inc_A;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Scalr_check( uplo, alpha, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  if ( FLA_Obj_equals( alpha, FLA_ONE ) )
  {
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_col_stride( A );
  inc_A    = 1;

  void* A_mat = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    A_mat = A_hip;
  }

  if ( uplo == FLA_LOWER_TRIANGULAR ){

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_alpha = ( float* ) FLA_FLOAT_PTR( alpha );
    float* buff_A_hip = ( float* ) A_mat;

    for ( i = 0; i < min( n_A, m_A ); i++ )
      rocblas_sscal( handle,
                     m_A - i,
                     buff_alpha,
                     buff_A_hip + i * ldim_A + i, inc_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_alpha = ( double* ) FLA_DOUBLE_PTR( alpha );
    double* buff_A_hip = ( double* ) A_mat;

    for ( i = 0; i < min( n_A, m_A ); i++ )
      rocblas_dscal( handle,
                     m_A - i,
                     buff_alpha,
                     buff_A_hip + i * ldim_A + i, inc_A );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_alpha = ( rocblas_float_complex* ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex* buff_A_hip = ( rocblas_float_complex* ) A_mat;

    for ( i = 0; i < min( n_A, m_A ); i++ )
      rocblas_cscal( handle,
                     m_A - i,
                     buff_alpha,
                     buff_A_hip + i * ldim_A + i, inc_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_alpha = ( rocblas_double_complex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex* buff_A_hip = ( rocblas_double_complex* ) A_mat;

    for ( i = 0; i < min( n_A, m_A ); i++ )
      rocblas_zscal( handle,
                     m_A - i,
                     buff_alpha,
                     buff_A_hip + i * ldim_A + i, inc_A );

    break;
  }

  }

  }

  else if ( uplo == FLA_UPPER_TRIANGULAR ){

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_alpha = ( float* ) FLA_FLOAT_PTR( alpha );
    float* buff_A_hip = ( float* ) A_mat;

    for ( i = 0; i < n_A; i++ )
      rocblas_sscal( handle,
                     min( i + 1, m_A ),
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_alpha = ( double* ) FLA_DOUBLE_PTR( alpha );
    double* buff_A_hip = ( double* ) A_mat;

    for ( i = 0; i < n_A; i++ )
      rocblas_dscal( handle,
                     min( i + 1, m_A ),
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_alpha = ( rocblas_float_complex* ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex* buff_A_hip = ( rocblas_float_complex* ) A_mat;

    for ( i = 0; i < n_A; i++ )
      rocblas_cscal( handle,
                     min( i + 1, m_A ),
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_alpha = ( rocblas_double_complex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex* buff_A_hip = ( rocblas_double_complex* ) A_mat;

    for ( i = 0; i < n_A; i++ )
      rocblas_zscal( handle,
                     min( i + 1, m_A ),
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A );

    break;
  }

  }

  }

  return FLA_SUCCESS;
}

#endif
