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

FLA_Error FLA_Trsm_external_hip( rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          ldim_A;
  int          ldim_B;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsm_check( side, uplo, trans, diag, alpha, A, B );

  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_length( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_length( B );

  rocblas_operation blas_trans = FLA_Param_map_flame_to_rocblas_trans( trans );
  rocblas_side blas_side = FLA_Param_map_flame_to_rocblas_side( side );
  rocblas_fill blas_uplo = FLA_Param_map_flame_to_rocblas_uplo( uplo );
  rocblas_diagonal blas_diag = FLA_Param_map_flame_to_rocblas_diag( diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    rocblas_strsm( handle,
                   blas_side,
                   blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( float * ) A_hip, ldim_A,
                   ( float * ) B_hip, ldim_B );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    rocblas_dtrsm( handle,
                 blas_side,
                 blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_B,
                 n_B,
                 buff_alpha,
                 ( double * ) A_hip, ldim_A,
                 ( double * ) B_hip, ldim_B );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_alpha = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( alpha );

    rocblas_ctrsm( handle,
                   blas_side,
                   blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( rocblas_float_complex * ) A_hip, ldim_A,
                   ( rocblas_float_complex * ) B_hip, ldim_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_alpha = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    rocblas_ztrsm( handle,
                   blas_side,
                   blas_uplo,
                   blas_trans,
                   blas_diag,
                   m_B,
                   n_B,
                   buff_alpha,
                   ( rocblas_double_complex * ) A_hip, ldim_A,
                   ( rocblas_double_complex * ) B_hip, ldim_B );

    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
