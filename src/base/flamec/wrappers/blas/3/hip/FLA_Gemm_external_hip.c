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

FLA_Error FLA_Gemm_external_hip( rocblas_handle handle, FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj beta, FLA_Obj C, void* C_hip )
{
  FLA_Datatype datatype;
  int          k_AB;
  int          m_A, n_A;
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
  ldim_A   = FLA_Obj_length( A );

  ldim_B   = FLA_Obj_length( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  ldim_C   = FLA_Obj_length( C );

  if ( transa == FLA_NO_TRANSPOSE || transa == FLA_CONJ_NO_TRANSPOSE )
    k_AB = n_A;
  else
    k_AB = m_A;

  rocblas_operation blas_transa = FLA_Param_map_flame_to_rocblas_trans( transa, FLA_Obj_is_real( A ) );
  rocblas_operation blas_transb = FLA_Param_map_flame_to_rocblas_trans( transb, FLA_Obj_is_real( B ) );


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
                   ( float * ) A_hip, ldim_A,
                   ( float * ) B_hip, ldim_B,
                   buff_beta,
                   ( float * ) C_hip, ldim_C );
    
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
                   ( double * ) A_hip, ldim_A,
                   ( double * ) B_hip, ldim_B,
                   buff_beta,
                   ( double * ) C_hip, ldim_C );
    
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
                   ( rocblas_float_complex * ) A_hip, ldim_A,
                   ( rocblas_float_complex * ) B_hip, ldim_B,
                   buff_beta,
                   ( rocblas_float_complex * ) C_hip, ldim_C );
    
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
                   ( rocblas_double_complex * ) A_hip, ldim_A,
                   ( rocblas_double_complex * ) B_hip, ldim_B,
                   buff_beta,
                   ( rocblas_double_complex * ) C_hip, ldim_C );
    
    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
