
#include "FLAME.h"

FLA_Error FLA_Chol_unb_external( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error    r_val = FLA_SUCCESS;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  int          info;
  FLA_Datatype datatype;
  int          m_A, cs_A;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Chol_check( uplo, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    F77_spotf2( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    F77_dpotf2( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    F77_cpotf2( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    F77_zpotf2( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  }

  // Convert to zero-based indexing, if an index was reported.
  if ( info > 0 ) r_val = info - 1;
  else            r_val = FLA_SUCCESS;

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return r_val;
}

FLA_Error FLA_Chol_l_unb_ext( FLA_Obj A )
{
  return FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, A );
}

FLA_Error FLA_Chol_u_unb_ext( FLA_Obj A )
{
  return FLA_Chol_unb_external( FLA_UPPER_TRIANGULAR, A );
}

