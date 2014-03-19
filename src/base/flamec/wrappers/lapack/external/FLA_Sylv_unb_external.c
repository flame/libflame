/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_unb_external( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          cs_A;
  int          cs_B;
  int          m_C, n_C, cs_C;
  char         blas_transa; 
  char         blas_transb;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Sylv_check( transa, transb, isgn, A, B, C, scale );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;
  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;
  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  cs_C     = FLA_Obj_col_stride( C );

  cs_A     = FLA_Obj_col_stride( A );

  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_netlib_trans( transa, &blas_transa );
  FLA_Param_map_flame_to_netlib_trans( transb, &blas_transb );


  switch( datatype ){

  case FLA_FLOAT:
  {
    int   *buff_isgn  = ( int   * ) FLA_INT_PTR( isgn );
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_C     = ( float * ) FLA_FLOAT_PTR( C );
    float *buff_scale = ( float * ) FLA_FLOAT_PTR( scale );

    F77_strsyl( &blas_transa,
                &blas_transb,
                buff_isgn,
                &m_C,
                &n_C,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_C, &cs_C,
                buff_scale,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    int    *buff_isgn  = ( int    * ) FLA_INT_PTR( isgn );
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_C     = ( double * ) FLA_DOUBLE_PTR( C );
    double *buff_scale = ( double * ) FLA_DOUBLE_PTR( scale );

    F77_dtrsyl( &blas_transa,
                &blas_transb,
                buff_isgn,
                &m_C,
                &n_C,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_C, &cs_C,
                buff_scale,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    int      *buff_isgn  = ( int      * ) FLA_INT_PTR( isgn );
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B     = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_C     = ( scomplex * ) FLA_COMPLEX_PTR( C );
    float    *buff_scale = ( float    * ) FLA_COMPLEX_PTR( scale );

    F77_ctrsyl( &blas_transa,
                &blas_transb,
                buff_isgn,
                &m_C,
                &n_C,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_C, &cs_C,
                buff_scale,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    int      *buff_isgn  = ( int      * ) FLA_INT_PTR( isgn );
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_C     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( C );
    double   *buff_scale = ( double   * ) FLA_DOUBLE_COMPLEX_PTR( scale );

    F77_ztrsyl( &blas_transa,
                &blas_transb,
                buff_isgn,
                &m_C,
                &n_C,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_C, &cs_C,
                buff_scale,
                &info );

    break;
  } 

  }

  // We don't provide a comprehensive strategy for handing scaling to avoid
  // overflow, so we just force the scale argument to 1.0.
  FLA_Set( FLA_ONE, scale );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error FLA_Sylv_nn_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_nh_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_hn_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
}

FLA_Error FLA_Sylv_hh_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
}

