/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Trsv_external_gpu( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu ) 
{
  FLA_Datatype datatype;
  integer          m_A;
  integer          ldim_A;
  integer          inc_x;
  char         blas_uplo;
  char         blas_trans;
  char         blas_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  ldim_A   = FLA_Obj_length( A );

  inc_x    = 1;

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );
  FLA_Param_map_flame_to_netlib_diag( diag, &blas_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    cublasStrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) x_gpu, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    cublasDtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) x_gpu, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    cublasCtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) x_gpu, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cublasZtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) x_gpu, inc_x );

    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
