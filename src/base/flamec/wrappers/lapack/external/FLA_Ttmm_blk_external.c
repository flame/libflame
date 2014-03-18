
#include "FLAME.h"

FLA_Error FLA_Ttmm_blk_external( FLA_Uplo uplo, FLA_Obj A )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, cs_A;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Ttmm_check( uplo, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    F77_slauum( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    F77_dlauum( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    F77_clauum( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    F77_zlauum( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error FLA_Ttmm_l_blk_ext( FLA_Obj A )
{
  return FLA_Ttmm_blk_external( FLA_LOWER_TRIANGULAR, A );
}

FLA_Error FLA_Ttmm_u_blk_ext( FLA_Obj A )
{
  return FLA_Ttmm_blk_external( FLA_UPPER_TRIANGULAR, A );
}

