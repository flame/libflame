/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Symm_external_gpu( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu )
{
  FLA_Datatype datatype;
  integer          m_C, n_C;
  integer          ldim_A;
  integer          ldim_B;
  integer          ldim_C;
  char         blas_side;
  char         blas_uplo; 

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Symm_check( side, uplo, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_length( A );

  ldim_B   = FLA_Obj_length( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  ldim_C   = FLA_Obj_length( C );

  FLA_Param_map_flame_to_netlib_side( side, &blas_side );
  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    cublasSsymm( blas_side,
                 blas_uplo,
                 m_C,
                 n_C,
                 *buff_alpha,
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) B_gpu, ldim_B,
                 *buff_beta,
                 ( float * ) C_gpu, ldim_C );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    cublasDsymm( blas_side,
                 blas_uplo,
                 m_C,
                 n_C,
                 *buff_alpha,
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) B_gpu, ldim_B,
                 *buff_beta,
                 ( double * ) C_gpu, ldim_C );
    
    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex *buff_alpha = ( cuComplex * ) FLA_COMPLEX_PTR( alpha );
    cuComplex *buff_beta  = ( cuComplex * ) FLA_COMPLEX_PTR( beta );

    cublasCsymm( blas_side,
                 blas_uplo,
                 m_C,
                 n_C,
                 *buff_alpha,
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) B_gpu, ldim_B,
                 *buff_beta,
                 ( cuComplex * ) C_gpu, ldim_C );
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex *buff_alpha = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    cuDoubleComplex *buff_beta  = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    cublasZsymm( blas_side,
                 blas_uplo,
                 m_C,
                 n_C,
                 *buff_alpha,
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) B_gpu, ldim_B,
                 *buff_beta,
                 ( cuDoubleComplex * ) C_gpu, ldim_C );
    
    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
