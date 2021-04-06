/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_unb_external( FLA_Obj A, FLA_Obj tu, FLA_Obj tv )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  integer          m_A, n_A, cs_A;
  integer          min_m_n, max_m_n;
  integer          lwork;
  FLA_Obj      d, e, work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Bidiag_check( A, tu, tv );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  min_m_n  = FLA_Obj_min_dim( A );
  max_m_n  = FLA_Obj_max_dim( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), min_m_n,     1, 0, 0, &d );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), min_m_n - 1, 1, 0, 0, &e );

  lwork    = max_m_n;
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float* buff_d    = ( float * ) FLA_FLOAT_PTR( d );
    float* buff_e    = ( float * ) FLA_FLOAT_PTR( e );
    float* buff_tu   = ( float * ) FLA_FLOAT_PTR( tu );
    float* buff_tv   = ( float * ) FLA_FLOAT_PTR( tv );
    float* buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );

    F77_sgebd2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_tu,
                buff_tv,
                buff_work,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double* buff_d    = ( double * ) FLA_DOUBLE_PTR( d );
    double* buff_e    = ( double * ) FLA_DOUBLE_PTR( e );
    double* buff_tu   = ( double * ) FLA_DOUBLE_PTR( tu );
    double* buff_tv   = ( double * ) FLA_DOUBLE_PTR( tv );
    double* buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

    F77_dgebd2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_tu,
                buff_tv,
                buff_work,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex* buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    float*    buff_d    = ( float    * ) FLA_FLOAT_PTR( d );
    float*    buff_e    = ( float    * ) FLA_FLOAT_PTR( e );
    scomplex* buff_tu   = ( scomplex * ) FLA_COMPLEX_PTR( tu );
    scomplex* buff_tv   = ( scomplex * ) FLA_COMPLEX_PTR( tv );
    scomplex* buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work_obj );

    F77_cgebd2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_tu,
                buff_tv,
                buff_work,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double*   buff_d    = ( double   * ) FLA_DOUBLE_PTR( d );
    double*   buff_e    = ( double   * ) FLA_DOUBLE_PTR( e );
    dcomplex* buff_tu   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( tu );
    dcomplex* buff_tv   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( tv );
    dcomplex* buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

    F77_zgebd2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_tu,
                buff_tv,
                buff_work,
                &info );

    break;
  } 

  }

  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error FLA_Bidiag_unb_ext( FLA_Obj A, FLA_Obj tu, FLA_Obj tv )
{
  return FLA_Bidiag_unb_external( A, tu, tv );
}
