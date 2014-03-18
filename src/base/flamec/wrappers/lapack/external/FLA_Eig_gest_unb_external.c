
#include "FLAME.h"

FLA_Error FLA_Eig_gest_unb_external( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Error    r_val = FLA_SUCCESS;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  int          itype;
  int          info;
  FLA_Datatype datatype;
  int          m_A, cs_A;
  int          cs_B;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

//  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  if ( inv == FLA_INVERSE )
    itype  = 1;
  else
    itype  = 2;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

    F77_ssygs2( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

    F77_dsygs2( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );

    F77_chegs2( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    F77_zhegs2( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  }

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return r_val;
}

FLA_Error FLA_Eig_gest_il_unb_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_unb_external( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_iu_unb_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_unb_external( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_nl_unb_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_unb_external( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_nu_unb_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_unb_external( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
}

