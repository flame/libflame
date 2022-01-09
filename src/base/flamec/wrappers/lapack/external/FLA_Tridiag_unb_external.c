/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_unb_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  integer          n_A, cs_A;
  FLA_Obj      d, e;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Tridiag_check( uplo, A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), n_A,     1, 0, 0, &d );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), n_A - 1, 1, 0, 0, &e );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float* buff_d    = ( float * ) FLA_FLOAT_PTR( d );
    float* buff_e    = ( float * ) FLA_FLOAT_PTR( e );
    float* buff_t    = ( float * ) FLA_FLOAT_PTR( t );

    F77_ssytd2( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double* buff_d    = ( double * ) FLA_DOUBLE_PTR( d );
    double* buff_e    = ( double * ) FLA_DOUBLE_PTR( e );
    double* buff_t    = ( double * ) FLA_DOUBLE_PTR( t );

    F77_dsytd2( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex* buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    float*    buff_d    = ( float    * ) FLA_FLOAT_PTR( d );
    float*    buff_e    = ( float    * ) FLA_FLOAT_PTR( e );
    scomplex* buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );

    F77_chetd2( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double*   buff_d    = ( double   * ) FLA_DOUBLE_PTR( d );
    double*   buff_e    = ( double   * ) FLA_DOUBLE_PTR( e );
    dcomplex* buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );

    F77_zhetd2( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                &info );

    break;
  } 

  }

  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error FLA_Tridiag_unb_ext( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  return FLA_Tridiag_unb_external( uplo, A, t );
}
