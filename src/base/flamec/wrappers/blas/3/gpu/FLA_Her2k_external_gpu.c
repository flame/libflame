
#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Her2k_external_gpu( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj B, void* B_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu )
{
  FLA_Datatype datatype;
  int          k_AB;
  int          m_A, n_A;
  int          m_C;
  int          ldim_A;
  int          ldim_B;
  int          ldim_C;
  char         blas_uplo; 
  char         blas_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Her2k_check( uplo, trans, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_length( A );

  ldim_B   = FLA_Obj_length( B );

  m_C      = FLA_Obj_length( C );
  ldim_C   = FLA_Obj_length( C );

  if ( trans == FLA_NO_TRANSPOSE )
    k_AB = n_A;
  else
    k_AB = m_A;

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );
    
    cublasSsyr2k( blas_uplo,
                  blas_trans,
                  m_C,
                  k_AB,
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
    
    cublasDsyr2k( blas_uplo,
                  blas_trans,
                  m_C,
                  k_AB,
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
    float *buff_beta = ( float * ) FLA_FLOAT_PTR( beta );
    
    cublasCher2k( blas_uplo,
                  blas_trans,
                  m_C,
                  k_AB,
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
    double *buff_beta = ( double * ) FLA_DOUBLE_PTR( beta );
    
    cublasZher2k( blas_uplo,
                  blas_trans,
                  m_C,
                  k_AB,
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
