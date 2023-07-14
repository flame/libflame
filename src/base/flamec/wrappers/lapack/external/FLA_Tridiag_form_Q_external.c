/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_form_Q_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  integer      info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  integer          m_A;
  integer          cs_A;
  integer          lwork;
  char         blas_uplo;
  FLA_Obj      work;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Tridiag_form_Q_check( uplo, A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

  lwork  = fla_max( 1, ( m_A - 1 ) ) * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work );
  

  switch( datatype ){

  case FLA_FLOAT:
  {
    float*    buff_A    = ( float    * ) FLA_FLOAT_PTR( A );
    float*    buff_t    = ( float    * ) FLA_FLOAT_PTR( t );
    float*    buff_work = ( float    * ) FLA_FLOAT_PTR( work );

    F77_sorgtr( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_t,
                buff_work, &lwork,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double*   buff_A    = ( double   * ) FLA_DOUBLE_PTR( A );
    double*   buff_t    = ( double   * ) FLA_DOUBLE_PTR( t );
    double*   buff_work = ( double   * ) FLA_DOUBLE_PTR( work );

    F77_dorgtr( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_t,
                buff_work, &lwork,
                &info );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex* buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex* buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );
    scomplex* buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work );

    F77_cungtr( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_t,
                buff_work, &lwork,
                &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
    dcomplex *buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work );

    F77_zungtr( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_t,
                buff_work, &lwork,
                &info );

    break;
  }

  }

  FLA_Obj_free( &work );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

