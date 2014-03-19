/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Gemv_external_gpu( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          ldim_A;
  int          inc_x;
  int          inc_y;
  char         blas_transa;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemv_check( transa, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_length( A );

  inc_x    = 1;
  inc_y    = 1;

  FLA_Param_map_flame_to_netlib_trans( transa, &blas_transa );


  switch( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    cublasSgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) x_gpu, inc_x,
                 *buff_beta,
                 ( float * ) y_gpu, inc_y );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    cublasDgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) x_gpu, inc_x,
                 *buff_beta,
                 ( double * ) y_gpu, inc_y );

    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex *buff_alpha = ( cuComplex * ) FLA_COMPLEX_PTR( alpha );
    cuComplex *buff_beta  = ( cuComplex * ) FLA_COMPLEX_PTR( beta );

    cublasCgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) x_gpu, inc_x,
                 *buff_beta,
                 ( cuComplex * ) y_gpu, inc_y );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex *buff_alpha = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    cuDoubleComplex *buff_beta  = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    cublasZgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) x_gpu, inc_x,
                 *buff_beta,
                 ( cuDoubleComplex * ) y_gpu, inc_y );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
