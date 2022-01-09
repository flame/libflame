/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Trmm_external_gpu( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu )
{
  FLA_Datatype datatype;
  integer          m_B, n_B;
  integer          ldim_A;
  integer          ldim_B;
  char         blas_side; 
  char         blas_uplo;
  char         blas_trans;
  char         blas_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trmm_check( side, uplo, trans, diag, alpha, A, B );

  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_length( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_length( B );

  FLA_Param_map_flame_to_netlib_side( side, &blas_side );
  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );
  FLA_Param_map_flame_to_netlib_diag( diag, &blas_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    cublasStrmm( blas_side,
                 blas_uplo, 
                 blas_trans,
                 blas_diag,
                 m_B,
                 n_B,
                 *buff_alpha,
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) B_gpu, ldim_B );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    cublasDtrmm( blas_side,
                 blas_uplo, 
                 blas_trans,
                 blas_diag,
                 m_B,
                 n_B,
                 *buff_alpha,
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) B_gpu, ldim_B );

    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex *buff_alpha = ( cuComplex * ) FLA_COMPLEX_PTR( alpha );

    cublasCtrmm( blas_side,
                 blas_uplo, 
                 blas_trans,
                 blas_diag,
                 m_B,
                 n_B,
                 *buff_alpha,
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) B_gpu, ldim_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex *buff_alpha = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    cublasZtrmm( blas_side,
                 blas_uplo, 
                 blas_trans,
                 blas_diag,
                 m_B,
                 n_B,
                 *buff_alpha,
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) B_gpu, ldim_B );

    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
