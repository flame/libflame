
#include "FLAME.h"

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Syrk_external_gpu( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj beta, FLA_Obj C, void* C_gpu )
{
  FLA_Datatype datatype;
  int          k_A;
  int          m_A, n_A;
  int          m_C;
  int          ldim_A;
  int          ldim_C;
  char         blas_uplo; 
  char         blas_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Syrk_check( uplo, trans, alpha, A, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_length( A );

  m_C      = FLA_Obj_length( C );
  ldim_C   = FLA_Obj_length( C );

  if ( trans == FLA_NO_TRANSPOSE )
    k_A = n_A;
  else
    k_A = m_A;

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    cublasSsyrk( blas_uplo,
                 blas_trans,
                 m_C,
                 k_A,
                 *buff_alpha,
                 ( float * ) A_gpu, ldim_A,
                 *buff_beta,
                 ( float * ) C_gpu, ldim_C );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    cublasDsyrk( blas_uplo,
                 blas_trans,
                 m_C,
                 k_A,
                 *buff_alpha,
                 ( double * ) A_gpu, ldim_A,
                 *buff_beta,
                 ( double * ) C_gpu, ldim_C );

    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex *buff_alpha = ( cuComplex * ) FLA_COMPLEX_PTR( alpha );
    cuComplex *buff_beta  = ( cuComplex * ) FLA_COMPLEX_PTR( beta );

    cublasCsyrk( blas_uplo,
                 blas_trans,
                 m_C,
                 k_A,
                 *buff_alpha,
                 ( cuComplex * ) A_gpu, ldim_A,
                 *buff_beta,
                 ( cuComplex * ) C_gpu, ldim_C );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex *buff_alpha = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    cuDoubleComplex *buff_beta  = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    cublasZsyrk( blas_uplo,
                 blas_trans,
                 m_C,
                 k_A,
                 *buff_alpha,
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 *buff_beta,
                 ( cuDoubleComplex * ) C_gpu, ldim_C );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
