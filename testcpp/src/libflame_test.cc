    /******************************************************************************
* Copyright (c) 2019 - present Advanced Micro Devices, Inc. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*******************************************************************************/

/*! @file liblame_test.cc
 *  libflame_test.cc Test application to validate CPP template interfaces
 *  for all libflame modules
 *  */

#include "libflame_test.hh"
typedef int FLA_Error;

FLA_Error geqp3_C( FLA_Obj A, int *jpvt, FLA_Obj t, FLA_Obj rwork )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;
  int          lwork;
  FLA_Obj      work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_QR_check( A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  lwork    = 3*n_A+1 ;//n_A * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
    sgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

   dgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  }

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
    lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
    lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
    float *r_work = ( float * ) FLA_FLOAT_PTR( rwork );

    cgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                r_work,
                &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
    lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
    lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
    double *r_work = ( double * ) FLA_DOUBLE_PTR( rwork );

    zgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                r_work,
                &info );

    break;
  }

  }

  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error potri_C( FLA_Uplo uplo, FLA_Obj A )
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

    spotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    dpotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );

    cpotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );

    zpotri_( &blas_uplo,
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

FLA_Error gelsd_C( FLA_Obj A, FLA_Obj B, FLA_Obj sCOBuff, FLA_Obj rcondC, int *rank )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A, cs_B;
  FLA_Obj work_obj, iwork_obj;
  FLA_Obj iwork_obj_float, iwork_obj_double;
  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_B     = FLA_Obj_col_stride( B );
  int NLVL = (m_A + n_A) *30;
  int lwork =12*m_A + 2*m_A*25 + 8*m_A*NLVL + m_A*n_A + (25+1)*2;

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );
  FLA_Obj_create( FLA_INT, lwork, 1, 0, 0, &iwork_obj );
  FLA_Obj_create( FLA_FLOAT, lwork, 1, 0, 0, &iwork_obj_float );
  FLA_Obj_create( FLA_DOUBLE, lwork, 1, 0, 0, &iwork_obj_double );
  switch( datatype ){

    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B    = ( float * ) FLA_FLOAT_PTR( B );
      float *buff_s    = ( float * ) FLA_FLOAT_PTR( sCOBuff );
      float *rcondCVal    =  ( float * ) FLA_FLOAT_PTR( rcondC );
      float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
      int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );

      sgelsd_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork, buff_iwork,
                  &info);

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B    = ( double * ) FLA_DOUBLE_PTR( B );
      double *buff_s    = ( double * ) FLA_DOUBLE_PTR( sCOBuff );
      double *rcondCVal  = ( double * ) FLA_DOUBLE_PTR( rcondC );
        double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );
      int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );

        dgelsd_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork, buff_iwork,
                  &info);

      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_B    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( B );
      float *buff_s    = ( float * ) FLA_FLOAT_PTR( sCOBuff );
      float *rcondCVal    =  ( float * ) FLA_FLOAT_PTR( rcondC );
      lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
      int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );
     float *buff_r    = ( float * ) FLA_FLOAT_PTR( iwork_obj_float );

     cgelsd_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal, rank,
                  buff_work, &lwork, buff_r,
                  buff_iwork, &info);

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_B    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( B );
      double *buff_s    = ( double * ) FLA_DOUBLE_PTR( sCOBuff );
      double *rcondCVal  = ( double * ) FLA_DOUBLE_PTR( rcondC );
        lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
      int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );
     double *buff_r    = ( double * ) FLA_DOUBLE_PTR( iwork_obj_double );
     zgelsd_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork, buff_r,
                  buff_iwork, &info);

      break;
    }
  }


  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error gelss_C( FLA_Obj A, FLA_Obj B, FLA_Obj sCOBuff, FLA_Obj rcondC, int *rank )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A, cs_B;
  FLA_Obj work_obj;
  FLA_Obj iwork_obj_float, iwork_obj_double;
  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_B     = FLA_Obj_col_stride( B );
  int NLVL = (m_A + n_A) *30;
  int lwork =12*m_A + 2*m_A*25 + 8*m_A*NLVL + m_A*n_A + (25+1)*2;

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );
  FLA_Obj_create( FLA_FLOAT, lwork, 1, 0, 0, &iwork_obj_float );
  FLA_Obj_create( FLA_DOUBLE, lwork, 1, 0, 0, &iwork_obj_double );
  switch( datatype ){

    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B    = ( float * ) FLA_FLOAT_PTR( B );
      float *buff_s    = ( float * ) FLA_FLOAT_PTR( sCOBuff );
      float *rcondCVal    =  ( float * ) FLA_FLOAT_PTR( rcondC );
      float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );

      sgelss_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork,
                  &info);

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B    = ( double * ) FLA_DOUBLE_PTR( B );
      double *buff_s    = ( double * ) FLA_DOUBLE_PTR( sCOBuff );
      double *rcondCVal  = ( double * ) FLA_DOUBLE_PTR( rcondC );
      double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

      dgelss_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork,
                  &info);

      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_B    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( B );
      float *buff_s    = ( float * ) FLA_FLOAT_PTR( sCOBuff );
      float *rcondCVal    =  ( float * ) FLA_FLOAT_PTR( rcondC );
      lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
      float *buff_r    = ( float * ) FLA_FLOAT_PTR( iwork_obj_float );

      cgelss_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal, rank,
                  buff_work, &lwork, buff_r,
                  &info);

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_B    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( B );
      double *buff_s    = ( double * ) FLA_DOUBLE_PTR( sCOBuff );
      double *rcondCVal  = ( double * ) FLA_DOUBLE_PTR( rcondC );
      lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
      double *buff_r    = ( double * ) FLA_DOUBLE_PTR( iwork_obj_double );
      zgelss_( &m_A,
                  &n_A,
                  &n_A,
                  buff_A, &cs_A,
                  buff_B, &cs_B,
                  buff_s, rcondCVal,
                  rank,
                  buff_work, &lwork, buff_r,
                  &info);

      break;
    }
  }


  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error larfg_C( int *n, FLA_Obj alphaObj, FLA_Obj A, int * incx, FLA_Obj tauObj)
{
  int          info = 0;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  datatype = FLA_Obj_datatype( A );

  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *alpha    = ( float * ) FLA_FLOAT_PTR( alphaObj );
      float *tau    = ( float * ) FLA_FLOAT_PTR( tauObj );
      slarfg_( n, alpha, buff_A, incx, tau);
      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *alpha    = ( double * ) FLA_DOUBLE_PTR( alphaObj );
      double *tau    = ( double * ) FLA_DOUBLE_PTR( tauObj );
      dlarfg_( n, alpha, buff_A, incx, tau);
      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *alpha    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( alphaObj );
      lapack_complex_float *tau    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( tauObj );
      clarfg_( n, alpha, buff_A, incx, tau);
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *alpha    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( alphaObj );
      lapack_complex_double *tau    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( tauObj );
      zlarfg_( n, alpha, buff_A, incx, tau);
      break;
    }
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

#if 0
FLA_Error larfgp_C( int *n, FLA_Obj alphaObj, FLA_Obj A, int * incx, FLA_Obj tauObj)
{
  int          info = 0;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  datatype = FLA_Obj_datatype( A );

  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *alpha    = ( float * ) FLA_FLOAT_PTR( alphaObj );
      float *tau    = ( float * ) FLA_FLOAT_PTR( tauObj );
      slarfgp_( n, alpha, buff_A, incx, tau);
      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *alpha    = ( double * ) FLA_DOUBLE_PTR( alphaObj );
      double *tau    = ( double * ) FLA_DOUBLE_PTR( tauObj );
      dlarfgp_( n, alpha, buff_A, incx, tau);
      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *alpha    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( alphaObj );
      lapack_complex_float *tau    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( tauObj );
      clarfgp_( n, alpha, buff_A, incx, tau);
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *alpha    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( alphaObj );
      lapack_complex_double *tau    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( tauObj );
      zlarfgp_( n, alpha, buff_A, incx, tau);
      break;
    }
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}
#endif

#if 0
FLA_Error orm2r_C( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  // int          m_A, n_A;
  int          m_B, n_B;
  int          cs_A;
  int          cs_B;
  int          k_t;
  int          lwork;
  char         blas_side;
  char         blas_trans;
  FLA_Obj      work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Apply_Q_check( side, trans, storev, A, t, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  // m_A      = FLA_Obj_length( A );
  // n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  cs_B     = FLA_Obj_col_stride( B );

  k_t      = FLA_Obj_vector_dim( t );

  FLA_Param_map_flame_to_netlib_side( side, &blas_side );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );

  if ( side == FLA_LEFT )
    lwork  = n_B * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  else
    lwork  = m_B * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float *buff_B    = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );

    if ( storev == FLA_COLUMNWISE )
      sorm2r_( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );
    else // storev == FLA_ROWWISE
      F77_sormlq( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double *buff_B    = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

    if ( storev == FLA_COLUMNWISE )
      F77_dorm2r( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );
    else // storev == FLA_ROWWISE
      F77_dormlq( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );
    scomplex *buff_B    = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work_obj );

    if ( storev == FLA_COLUMNWISE )
      F77_cunm2r( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );
    else // storev == FLA_ROWWISE
      F77_cunmlq( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
    dcomplex *buff_B    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

    if ( storev == FLA_COLUMNWISE )
      F77_zunm2r( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );
    else // storev == FLA_ROWWISE
      F77_zunmlq( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );

    break;
  }

  }

  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}
#endif

FLA_Error FLA_Sylv_unb_external1( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
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

    strsyl_( &blas_transa,
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

    dtrsyl_( &blas_transa,
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

 // case FLA_COMPLEX:
 // {
 //   int      *buff_isgn  = ( int      * ) FLA_INT_PTR( isgn );
 //   lapack_complex_float *buff_A     = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
 //   lapack_complex_float *buff_B     = ( lapack_complex_float * ) FLA_COMPLEX_PTR( B );
 //   lapack_complex_float *buff_C     = ( lapack_complex_float * ) FLA_COMPLEX_PTR( C );
 //   float    *buff_scale = ( float    * ) FLA_COMPLEX_PTR( scale );
 //
 //   ctrsyl_( &blas_transa,
 //               &blas_transb,
 //               buff_isgn,
 //               &m_C,
 //               &n_C,
 //               buff_A, &cs_A,
 //               buff_B, &cs_B,
 //               buff_C, &cs_C,
 //               buff_scale,
 //               &info );
 //
 //   break;
 // }
 //
 // case FLA_DOUBLE_COMPLEX:
 // {
 //   int      *buff_isgn  = ( int      * ) FLA_INT_PTR( isgn );
 //   lapack_complex_double *buff_A     = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
 //   lapack_complex_double *buff_B     = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( B );
 //   lapack_complex_double *buff_C     = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( C );
 //   double   *buff_scale = ( double   * ) FLA_DOUBLE_COMPLEX_PTR( scale );
 //
 //   ztrsyl_( &blas_transa,
 //               &blas_transb,
 //               buff_isgn,
 //               &m_C,
 //               &n_C,
 //               buff_A, &cs_A,
 //               buff_B, &cs_B,
 //               buff_C, &cs_C,
 //               buff_scale,
 //               &info );
 //
 //   break;
 // }
 //
 }

  // We don't provide a comprehensive strategy for handing scaling to avoid
  // overflow, so we just force the scale argument to 1.0.
  //FLA_Set( FLA_ONE, scale );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

template< typename T >
void potrf_test()
{
  int m = 32 ;
  char uplo = 'U';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potrf(LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  //Call C function
  FLA_Chol_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj );
  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "potrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}


template< typename T >
void potf2_test()
{

  int m = 16; //2048
  char uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potf2( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potf2():Failure Diff = %E\n", diff);
  }else{
    printf( "potf2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void getrf_test()
{

  int m = 512;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, pivCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int *pivCPPOBuff, *pivCOBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  pivCPPOBuff =  new int [min_m_n];
  pivCOBuff =  new int [min_m_n];


 //Call CPP function
  libflame::getrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, pivCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( pivCOBuff, 1, min_m_n, &pivCOObj );
  //FLA_Set( FLA_ZERO, pivCOObj );

  //Call C function
  FLA_LU_piv_blk_external( aCIOObj, pivCOObj );


  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivCOBuff, pivCPPOBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getrf(): Failure Pivot BUffer mismatach Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete pivCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &pivCOObj );
}


template< typename T >
void getf2_test()
{

  int m = 1024;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, pivCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int *pivCPPOBuff, *pivCOBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  pivCPPOBuff =  new int [min_m_n];
  pivCOBuff =  new int [min_m_n];

 //Call CPP function
  libflame::getf2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, pivCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( pivCOBuff, 1, min_m_n, &pivCOObj );
  //FLA_Set( FLA_ZERO, pivCOObj );

  //Call C function
  FLA_LU_piv_unb_external( aCIOObj, pivCOObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivCOBuff, pivCPPOBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getf2(): Failure:Pivot Buffer Mismatch Diff = %d\n", diffInt);
  }else{
    printf( "getf2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete pivCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &pivCOObj );
}
#if 0
#define PREFIX2FLAME_INVERT_TAU(datatype, tauCOObj)\
{\
  if (datatype == FLA_FLOAT)\
    FLAME_invert_stau(tauCOObj);\
  if (datatype == FLA_DOUBLE)\
    FLAME_invert_dtau(tauCOObj);\
  if (datatype == FLA_COMPLEX)\
    FLAME_invert_ctau(tauCOObj);\
  if (datatype == FLA_DOUBLE_COMPLEX)\
    FLAME_invert_ztau(tauCOObj);\
}
#endif

template< typename T >
void geqrf_test()
{

  int m = 512;
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int min_m_n = min( m, n );
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

 //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );

  double diff = computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "geqrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "geqrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void geqr2_test()
{


  int m = 128;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int min_m_n = min( m, n ) ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::geqr2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call C function
  FLA_QR_unb_external( aCIOObj, tauCOObj );



  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqr2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqr2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void geqpf_test()
{
  int m = 256;
  int n = 256;//4096;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  T *aCPPIOBuff, *aCIOBuff, *rworkRefBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  rworkRefBuff = new T [2 * n];

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(jpvtCPPOBuff, jpvtCOBuff, n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatype, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff );

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "geqpf(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}

template< typename Ta, typename Tb >
void geqpf_test()
{

  int m = 128;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new Ta [min_m_n];
  tauCOBuff =  new Ta [min_m_n];
  rworkRefBuff = new Tb [2 * n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(jpvtCPPOBuff, jpvtCOBuff, n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;
  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff ) ;

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );


  double diff = computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff += computeError<Ta>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  int diffInt = computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "geqpf(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}



template< typename T >
void geqp3_test()
{
  int m = 256;
  int n = 4096;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  T *aCPPIOBuff, *aCIOBuff, *rworkRefBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  rworkRefBuff = new T [2 * n];

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(jpvtCPPOBuff, jpvtCOBuff, n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatype, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff );

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );


  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "geqp3(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}

template< typename Ta, typename Tb >
void geqp3_test()
{
  int m = 128;//2048;
  int n = 128;//512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new Ta [min_m_n];
  tauCOBuff =  new Ta [min_m_n];
  rworkRefBuff = new Tb [2 * n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(jpvtCPPOBuff, jpvtCOBuff, n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff ) ;

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );

  double diff = computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff += computeError<Ta>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  int diffInt = computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "geqp3(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff1111 = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}


template< typename T >
void gelqf_test()
{
  int m = 8192;
  int n = 1024;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int min_m_n = min( m, n ) ;

  int datatype = getDatatype<T>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );



  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelqf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelqf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gelq2_test()
{
  int m = 1000;
  int n = 1000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gelq2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &lda, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call C function
  FLA_LQ_unb_external( aCIOObj, tauCOObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelq2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelq2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

extern TLS_CLASS_SPEC FLA_Obj FLA_ZERO;
template< typename T >
void gelsd_test()
{
  int m = 128;
  int n = 1024;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj, sCOObj;
  T *aCPPIOBuff, *bCPPIOBuff, *aCIOBuff, *bCIOBuff ;
  T  *sCPPOBuff, *sCOBuff  ;
  T rcondCPP = 0;
  int rankCPP ;
  int rankC ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));


  sCPPOBuff = new T [min_m_n];
  sCOBuff = new T [min_m_n];

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, ldb*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, ldb, n, &bCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, ldb, &bCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );

  //Call CPP function
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, &n, aCPPIOBuff, &lda, bCPPIOBuff, &ldb, sCPPOBuff, &rcondCPP, &rankCPP );

  //Call C function
  gelsd_C( aCIOObj, bCIOObj, sCOObj, FLA_ZERO, &rankC );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff +=  computeError<T>( ldb, n, bCIOBuff, bCPPIOBuff ) ;
  diff +=  computeError<T>( 1, min_m_n, sCOBuff, sCPPOBuff ) ;
  diff += abs(rankCPP - rankC);

  if(diff != 0.0)
  {
    printf( "gelsd(): Failure Diff = %E\n", diff) ;
  }else{
    printf( "gelsd(): Success\n");
  }

 //Free up the buffers
 delete aCPPIOBuff ;
 delete bCPPIOBuff ;
 delete sCPPOBuff ;
 FLA_Obj_free( &aCIOObj );
 FLA_Obj_free( &bCIOObj );
 FLA_Obj_free(&sCOObj );
}

template< typename Ta, typename Tb >
void gelsd_test()
{
  int m = 16384;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj, sCOObj;
  Ta *aCPPIOBuff, *bCPPIOBuff, *aCIOBuff, *bCIOBuff ;
  Tb  *sCPPOBuff, *sCOBuff  ;
  Tb rcondCPP = 0;
  int rankCPP ;
  int rankC ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));


  sCPPOBuff = new Tb [min_m_n];
  sCOBuff = new Tb [min_m_n];

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, ldb*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, ldb, n, &bCIOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, min_m_n, 1, &sCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, ldb, &bCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );

  //Call CPP function
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, &n, aCPPIOBuff, &lda, bCPPIOBuff, &ldb, sCPPOBuff, &rcondCPP, &rankCPP );

  //Call C function
  gelsd_C( aCIOObj, bCIOObj, sCOObj, FLA_ZERO, &rankC );

  double diff =  computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff +=  computeError<Ta>( ldb, n, bCIOBuff, bCPPIOBuff ) ;
  diff +=  computeError<Tb>( 1, min_m_n, sCOBuff, sCPPOBuff ) ;
  diff += abs(rankCPP - rankC);

  if(diff != 0.0)
  {
    printf( "gelsd(): Failure Diff = %E\n", diff) ;
  }else{
    printf( "gelsd(): Success\n");
  }

 //Free up the buffers
 delete aCPPIOBuff ;
 delete bCPPIOBuff ;
 delete sCPPOBuff ;
 FLA_Obj_free( &aCIOObj );
 FLA_Obj_free( &bCIOObj );
 FLA_Obj_free(&sCOObj );
}

extern TLS_CLASS_SPEC FLA_Obj FLA_ZERO;
template< typename T >
void gelss_test()
{
  int m = 128;
  int n = 1024;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj, sCOObj;
  T *aCPPIOBuff, *bCPPIOBuff, *aCIOBuff, *bCIOBuff ;
  T  *sCPPOBuff, *sCOBuff  ;
  T rcondCPP = 0;
  int rankCPP ;
  int rankC ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));


  sCPPOBuff = new T [min_m_n];
  sCOBuff = new T [min_m_n];

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, ldb*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, ldb, n, &bCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, ldb, &bCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );

  //Call CPP function
  libflame::gelss( LAPACK_COL_MAJOR, &m, &n, &n, aCPPIOBuff, &lda, bCPPIOBuff, &ldb, sCPPOBuff, &rcondCPP, &rankCPP );

  //Call C function
  gelss_C( aCIOObj, bCIOObj, sCOObj, FLA_ZERO, &rankC );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff +=  computeError<T>( ldb, n, bCIOBuff, bCPPIOBuff ) ;
  diff +=  computeError<T>( 1, min_m_n, sCOBuff, sCPPOBuff ) ;
  diff += abs(rankCPP - rankC);

  if(diff != 0.0)
  {
    printf( "gelss(): Failure Diff = %E\n", diff) ;
  }else{
    printf( "gelss(): Success\n");
  }

 //Free up the buffers
 delete aCPPIOBuff ;
 delete bCPPIOBuff ;
 delete sCPPOBuff ;
 FLA_Obj_free( &aCIOObj );
 FLA_Obj_free( &bCIOObj );
 FLA_Obj_free(&sCOObj );
}

template< typename Ta, typename Tb >
void gelss_test()
{
  int m = 16384;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj, sCOObj;
  Ta *aCPPIOBuff, *bCPPIOBuff, *aCIOBuff, *bCIOBuff ;
  Tb  *sCPPOBuff, *sCOBuff  ;
  Tb rcondCPP = 0;
  int rankCPP ;
  int rankC ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));


  sCPPOBuff = new Tb [min_m_n];
  sCOBuff = new Tb [min_m_n];

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, ldb*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, ldb, n, &bCIOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, min_m_n, 1, &sCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, ldb, &bCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );

  //Call CPP function
  libflame::gelss( LAPACK_COL_MAJOR, &m, &n, &n, aCPPIOBuff, &lda, bCPPIOBuff, &ldb, sCPPOBuff, &rcondCPP, &rankCPP );

  //Call C function
  gelss_C( aCIOObj, bCIOObj, sCOObj, FLA_ZERO, &rankC );

  double diff =  computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff +=  computeError<Ta>( ldb, n, bCIOBuff, bCPPIOBuff ) ;
  diff +=  computeError<Tb>( 1, min_m_n, sCOBuff, sCPPOBuff ) ;
  diff += abs(rankCPP - rankC);

  if(diff != 0.0)
  {
    printf( "gelss(): Failure Diff = %E\n", diff) ;
  }else{
    printf( "gelss(): Success\n");
  }

 //Free up the buffers
 delete aCPPIOBuff ;
 delete bCPPIOBuff ;
 delete sCPPOBuff ;
 FLA_Obj_free( &aCIOObj );
 FLA_Obj_free( &bCIOObj );
 FLA_Obj_free(&sCOObj );
}

template< typename T >
void lauum_test()
{

  int m = 64;
  char uplo = 'u';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::lauum( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  FLA_Ttmm_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "lauum(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauum(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void lauu2_test()
{

  int m = 128;
  char uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::lauu2( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Param_map_flame_to_netlib_uplo( FLA_LOWER_TRIANGULAR, &uplo );

  //Call C function
  FLA_Ttmm_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "lauu2(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauu2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}


template< typename T >
void potri_test()
{
  int m = 64; //64 works
  char uplo = 'L';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potri( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  potri_C( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potri(): Failure Diff = %E\n", diff);
  }else{
    printf( "potri(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void trtri_test()
{
  int m = 128;
  char uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::trtri( LAPACK_COL_MAJOR, &uplo, &blas_diag, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trtri(): Failure Diff = %E\n", diff);
  }else{
    printf( "trtri(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void trti2_test()
{
  int m = 128;
  char uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::trti2( LAPACK_COL_MAJOR, &uplo, &blas_diag, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );
  //FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trti2(): Failure Diff = %E\n", diff);
  }else{
    printf( "trti2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}
extern TLS_CLASS_SPEC FLA_Obj FLA_ONE;
template< typename T >
void trsyl_test()
{
  int m = 64;
  int n = 64;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  FLA_Trans transARef = FLA_NO_TRANSPOSE ; FLA_Trans transBRef = FLA_NO_TRANSPOSE;
  FLA_Obj aCIOObj,  bCIOObj, cCIOObj, scaleRefObj ;
  T *aCPPIOBuff, *bCPPIOBuff, *cCPPIOBuff ;
  T *aCIOBuff, *bCIOBuff, *cCIOBuff, scaleRef ;
  T scale;
  int isgn = 1;
  int isgnValue = 1;
  int isgnRefValue  = 1;
  int *isgnRef = &isgnValue;

   FLA_Obj TLS_CLASS_SPEC isgnObj = {};

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);

//int matrix_layout, char* transa, char* transb, int* isgn, int* m, int* n, T* a, int* lda, T* b, int* ldb, T* c, int* ldc, T* scale )
  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bCPPIOBuff, &n, cCPPIOBuff, &m, &scale );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, n, &bCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( &scaleRef, 1, 1, &scaleRefObj );

  //Call C function
  FLA_Sylv_unb_external( transARef, transBRef, FLA_ONE, aCIOObj,  bCIOObj, cCIOObj, scaleRefObj );


  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): Failure Diff = %E\n", diff);
  }else{
   printf( "trsyl(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bCPPIOBuff;
  delete cCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
  FLA_Obj_free( &cCIOObj );
}

template< typename Ta, typename Tb  >
void trsyl_test()
{
  int m = 64;
  int n = 64;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  FLA_Trans transARef = FLA_NO_TRANSPOSE ; FLA_Trans transBRef = FLA_NO_TRANSPOSE;
  FLA_Obj aCIOObj,  bCIOObj, cCIOObj, scaleRefObj ;
  Tb *aCPPIOBuff, *bCPPIOBuff, *cCPPIOBuff ;
  Tb *aCIOBuff, *bCIOBuff, *cCIOBuff, scaleRef ;
  Tb scale;
  int isgn = 1;
  int isgnValue = 1;

  int *isgnRef = &isgnValue;
   FLA_Obj TLS_CLASS_SPEC isgnObj = {};


  int datatype = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);

  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bCPPIOBuff, &n, cCPPIOBuff, &m, &scale );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, n, &bCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( &scaleRef, 1, 1, &scaleRefObj );

  //Call C function
  FLA_Sylv_unb_external( transARef, transBRef, FLA_ONE, aCIOObj,  bCIOObj, cCIOObj, scaleRefObj );


  //Compute Difference in C and CPP buffer
  double diff =  computeError<Tb>( m, m, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): Failure Diff = %E\n", diff);
  }else{
   printf( "trsyl(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bCPPIOBuff;
  delete cCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
  FLA_Obj_free( &cCIOObj );
}

template< typename T >
void gehrd_test()
{
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  tauCOBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::gehrd( LAPACK_COL_MAJOR, &n, &ilo, &iho, aCPPIOBuff, &n, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Hess_blk_external( aCIOObj, tauCOObj, ilo, iho );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gehrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gehd2_test()
{
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);

  tauCPPOBuff =  new T [n-1];
  tauCOBuff =  new T [n-1];
  for(int j=0; j < n-1; j++)
  {
    tauCPPOBuff[j] = 0;
    tauCOBuff[j] = 0;
  }

  //Call CPP function
  libflame::gehd2( LAPACK_COL_MAJOR, &n, &ilo, &iho, aCPPIOBuff, &n, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Hess_unb_external( aCIOObj, tauCOObj, ilo, iho );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gehd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void sytrd_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *d, *e ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauCOBuff =  new T [n-1];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "sytrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta, typename Tb >
void hetrd_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb  *d, *e ;
  char uplo = 'u';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauCOBuff =  new Ta [n-1];

  //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "hetrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void sytd2_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff, *d, *e ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauCOBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::sytd2( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "sytd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta,  typename Tb >
void hetd2_test()
{
  int n = 128;//2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *d, *e ;
  char uplo = 'L';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauCOBuff =  new Ta [n-1];

 for (int j=0; j < n-1; j++)
 {
    tauCPPOBuff[j] = 0;
    tauCOBuff[j] = 0;
    e[j] = 0;
 }
  for (int j=0; j < n; j++)
 {
    d[j] = 0;
 }

  //Call CPP function
  libflame::hetd2( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "hetd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gebrd_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<T>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void gebrd_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void gebd2_test()
{

  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_unb_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<T>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}


template< typename Ta, typename Tb >
void gebd2_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];
   fp = fopen("test/in","a+");
   print(fp,m*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_unb_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupCOBuff, taupCPPOBuff );
   fp = fopen("test/tau1","a+");
   print(fp,min_m_n,tauqCOBuff, taupCOBuff);
   fclose(fp);
   fp = fopen("test/tau2","a+");
   print(fp,min_m_n,tauqCPPOBuff, taupCPPOBuff);
   fclose(fp);
     fp = fopen("test/out","a+");
   print(fp,m*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
     fp = fopen("test/out1","a+");
   print(fp,m*n,aCPPIOBuff, aCPPIOBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,m*n,aCIOBuff, aCIOBuff);
   fclose(fp);

  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void sygst_test()
{
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);


  //Call CPP function
  libflame::sygst( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bCIOObj );

  //Call C function
  FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bCIOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygst(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
}

template< typename T >
void hegst_test()
{

  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);


  //Call CPP function
  libflame::hegst( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bCIOObj );

  //Call C function
  FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bCIOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegst(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
}

template< typename T >
void sygs2_test()
{
  int n = 64;//2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);


  //Call CPP function
  libflame::sygs2( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bCIOObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bCIOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
}

template< typename T >
void hegs2_test()
{
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);


  //Call CPP function
  libflame::hegs2( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bCIOObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bCIOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegs2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
}

template< typename T >
void larft_test()
{
#if 0
  int n = 64;//256;
  int k = 128;//128;
 int min_m_n = min( m, n );
 srand (time(NULL));

 FLA_Init( );
 FLA_Obj vCIObj, tauCIObj, tCOObj;
 T *vCPPIBuff, *vCIBuff, *tCPPIBuff, *tCOBuff ;
 T *tauCPPIBuff, *tauCIBuff ;
 int datatype = getDatatype<T>();
 char direct  = 'F';
 char storev = 'C';
 int ldv = n;
 int ldt = k;
 //Allocate and initialize buffers for C and CPP functions with random values
 allocate_init_buffer(vCPPIBuff, vCIBuff, ldv*k);
 tauCPPIBuff =  new T [min_m_n];
 tauCIBuff =  new T [min_m_n];
 tCPPOBuff =  new T [ldt*k];
 tCOBuff =  new T [ldt*k];

 //Call CPP function
  libflame::larft( LAPACK_COL_MAJOR, &direct, &storev, &n, &k, &vCPPIBuff, &ldv, tauCPPIBuff, tCPPOBuff, *ldt );

//int larft( int matrix_layout, char* direct, char* storev, int* n, int* k, T* v, int* ldv, T* tau, T* t, int* ldt )
//{
//  return larft( matrix_layout, direct, storev, n, k, v, ldv, tau, t, ldt );
//}

 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, ldv, k, &vCIObj );
 FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCIObj );
 FLA_Obj_create_without_buffer( datatype, ldt, k, &tCOObj );
 FLA_Obj_attach_buffer( vCIBuff, 1, ldv, &vCIObj );
 FLA_Obj_attach_buffer( tauCIBuff, 1, min_m_n, &tauCIObj );
 FLA_Obj_attach_buffer( tCOBuff, 1, ldt, &tCOObj );

  //Call C function
  larft_C( aCIOObj, tauCOObj );
  double diff =  computeError<T>( ldt, k, tCOBuff, tCPPOBuff );

  if(diff != 0.0)
  {
    printf( "larft(): Failure Diff = %E\n", diff);
  }else{
    printf( "larft(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  #endif
}

template< typename T >
void larfg_test()
{

  int n = 256;
  srand (time(NULL));

  FLA_Init( );
  T alphaC = 1, alphaCPP =1;
  T tauC, tauCPP;
  FLA_Obj xCIOObj, alphaCIOObj, tauCOObj;
  T *xCPPIOBuff, *xCIOBuff;

  int datatype = getDatatype<T>();
  int incx = 5;
  int size  = (1+(n-2)*abs(incx));
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(xCPPIOBuff, xCIOBuff, size);

  //Call CPP function
  libflame::larfg( &n, &alphaCPP, xCPPIOBuff, &incx, &tauCPP);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, size, 1, &xCIOObj );
  FLA_Obj_attach_buffer( xCIOBuff, 1, size, &xCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &alphaCIOObj );
  FLA_Obj_attach_buffer( &alphaC, 1, 1, &alphaCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( &tauC, 1, 1, &tauCOObj );

  //Call C function
  larfg_C( &n, alphaCIOObj, xCIOObj, &incx, tauCOObj);
  double diff =  computeError<T>( size, 1, xCIOBuff, xCPPIOBuff );
  diff += computeError<T>( 1, 1, &tauC, &tauCPP);
  diff += computeError<T>( 1, 1, &alphaC, &alphaCPP);

  if(diff != 0.0)
  {
    printf( "larfg(): Failure Diff = %E\n", diff);
  }
  else{
    printf( "larfg(): Success\n");
  }

  //Free up the buffers
  delete xCPPIOBuff ;
  FLA_Obj_free( &xCIOObj );
}


template< typename T >
void larfgp_test()
{

  int n = 256;
  srand (time(NULL));

  FLA_Init( );
  T alphaC = 1, alphaCPP =1;
  T tauC, tauCPP;
  FLA_Obj xCIOObj, alphaCIOObj, tauCOObj;
  T *xCPPIOBuff, *xCIOBuff;

  int datatype = getDatatype<T>();
  int incx = 5;
  int size  = (1+(n-2)*abs(incx));
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(xCPPIOBuff, xCIOBuff, size);

  //Call CPP function
  libflame::larfgp( LAPACK_COL_MAJOR, &n, &alphaCPP, xCPPIOBuff, &incx, &tauCPP);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, size, 1, &xCIOObj );
  FLA_Obj_attach_buffer( xCIOBuff, 1, size, &xCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &alphaCIOObj );
  FLA_Obj_attach_buffer( &alphaC, 1, 1, &alphaCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( &tauC, 1, 1, &tauCOObj );

  //Call C function
 // larfgp_C( &n, alphaCIOObj, xCIOObj, &incx, tauCOObj);
  double diff =  computeError<T>( size, 1, xCIOBuff, xCPPIOBuff );
  diff += computeError<T>( 1, 1, &tauC, &tauCPP);
  diff += computeError<T>( 1, 1, &alphaC, &alphaCPP);

  if(diff != 0.0)
  {
    printf( "larfgp(): Failure Diff = %E\n", diff);
  }
  else{
    printf( "larfgp(): Success\n");
  }

  //Free up the buffers
  delete xCPPIOBuff ;
  FLA_Obj_free( &xCIOObj );
}


template< typename T >
void orgqr_test()
{
  int m = 64;//512;
  int n = 64;//256;
  int k = 64;//128;
 int min_m_n = min( m, n );
 srand (time(NULL));

 FLA_Init( );
 FLA_Obj aCIOObj, tauCOObj;
 T *aCPPIOBuff, *aCIOBuff ;
 T *tauCPPOBuff, *tauCOBuff ;
 int datatype = getDatatype<T>();

 //Allocate and initialize buffers for C and CPP functions with random values
 allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
 tauCPPOBuff =  new T [min_m_n];
 tauCOBuff =  new T [min_m_n];

 //Call CPP function
 libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
 libflame::orgqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
 FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
 FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
 FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_QR_form_Q_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orgqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ungqr_test()
{
  int m = 512;
  int n = 256;
  int k = 128;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
  libflame::ungqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_QR_form_Q_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "ungqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ormqr_test()
{
  int m = 64;//512;
  int n = 128;//256;
  int k = 64;//128;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'T';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*k);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, aCPPIOBuff, &n, tauCPPOBuff );
  libflame::ormqr( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &n, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, k, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void unmqr_test()
{
  int m = 512;
  int n = 1024;
  int k = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'N';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*k);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, aCPPIOBuff, &n, tauCPPOBuff );
  libflame::unmqr( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &n, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, k, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unmqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void orm2r_test()
{
  int m = 64;//512;
  int n = 128;//256;
  int k = 64;//128;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'T';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*k);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, aCPPIOBuff, &n, tauCPPOBuff );
  libflame::orm2r( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &n, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, k, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  //orm2r_C( FLA_RIGHT, FLA_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orm2r(): Failure Diff = %E\n", diff);
  }else{
    printf( "orm2r(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void unm2r_test()
{
  int m = 512;
  int n = 1024;
  int k = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'N';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*k);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, aCPPIOBuff, &n, tauCPPOBuff );
  libflame::unm2r( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &n, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, k, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
 // orm2r_C( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unm2r(): Failure Diff = %E\n", diff);
  }else{
    printf( "unm2r(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void orglq_test()
{
  int m = 512;
  int n = 256;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *work;
  char uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [k-1];
  tauCOBuff =  new T [k-1];
  work = new T [lwork];

  //Call CPP function
  libflame::orglq( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k-1, &tauCOObj );

  //Call C function
 // orglq_C( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orglq(): Failure Diff = %E\n", diff);
  }else{
    printf( "orglq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ormlq_test()
{
  int m = 512;
  int n = 256;
  int k = 128;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'T';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, k*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &k, &n, aCPPIOBuff, &k, tauCPPOBuff );
  libflame::ormlq( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &k, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, k, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, k, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_TRANSPOSE, FLA_ROWWISE, aCIOObj, tauCOObj, cCIOObj );

  double diff =  computeError<T>( k, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormlq(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormlq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void unmlq_test()
{
  int m = 256;
  int n = 512;
  int k = 128;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'N';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, k*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &k, &n, aCPPIOBuff, &k, tauCPPOBuff );
  libflame::unmlq( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &k, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, k, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, k, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_ROWWISE, aCIOObj, tauCOObj, cCIOObj );

  double diff =  computeError<T>( k, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unmlq(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmlq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void orgtr_test()
{
  int m = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  tauCPPOBuff =  new T [m-1];
  tauCOBuff =  new T [m-1];
  T *d =  new T [m];
  T *e =  new T [m-1];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::orgtr( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, tauCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_form_Q_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, m-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orgtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta,  typename Tb>
void ungtr_test()
{
  int m = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  char uplo = 'U';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  tauCPPOBuff =  new Ta [m-1];
  tauCOBuff =  new Ta [m-1];
  Tb *d =  new Tb [m];
  Tb *e =  new Tb [m-1];

  //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::ungtr( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, tauCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_form_Q_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( m, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, m-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "ungtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ormtr_test()
{
  int m = 512;
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  char sideCPP = 'L';
  char uplo = 'L';
  char transCPP = 'T';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [m-1];
  tauCOBuff =  new T [m-1];
  T *d =  new T [m];
  T *e =  new T [m-1];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::ormtr( LAPACK_COL_MAJOR, &sideCPP, &uplo, &transCPP, &m, &n, aCPPIOBuff, &m, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_apply_Q_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta, typename Tb >
void unmtr_test()
{
  int m = 8;
  int n = 8;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  char sideCPP = 'R';
  char uplo = 'U';
  char transCPP = 'N';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new Ta [m-1];
  tauCOBuff =  new Ta [m-1];
  Tb *d =  new Tb [m];
  Tb *e =  new Tb [m-1];

    //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::unmtr( LAPACK_COL_MAJOR, &sideCPP, &uplo, &transCPP, &m, &n, aCPPIOBuff, &m, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_apply_Q_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unmtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete cCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &cCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void orgbr_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  char vect = 'P';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::orgbr( LAPACK_COL_MAJOR, &vect, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_form_V_external( aCIOObj, taupCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orgbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void ungbr_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  char vect = 'P';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::ungbr( LAPACK_COL_MAJOR, &vect, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_form_V_external( aCIOObj, taupCOObj );

  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ungbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void ormbr_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  char vect = 'Q';
  char sideCPP = 'L';
  char transCPP = 'T';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::ormbr( LAPACK_COL_MAJOR, &vect, &sideCPP, &transCPP, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff, cCPPIOBuff, &m  );
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_apply_U_external( FLA_LEFT, FLA_TRANSPOSE, aCIOObj, taupCOObj, cCIOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void unmbr_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *cCPPIOBuff, *cCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  char vect = 'P';
  char sideCPP = 'R';
  char transCPP = 'N';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::unmbr( LAPACK_COL_MAJOR, &vect, &sideCPP, &transCPP, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff, cCPPIOBuff, &m  );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj );
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj );

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_apply_V_external( FLA_RIGHT, FLA_NO_TRANSPOSE, aCIOObj, taupCOObj, cCIOObj );

  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unmbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void steqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  T *dCPPIOBuff, *eCPPIOBuff, *zCPPIOBuff;
  T *dCIOBuff, *eCIOBuff, *zCIOBuff;
  char jobz = 'N';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::steqr( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevd_external( FLA_EVD_WITHOUT_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "steqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "steqr(): Success\n");
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );

}

template< typename Ta, typename Tb >
void steqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  Tb *dCPPIOBuff, *eCPPIOBuff;
  Ta *zCPPIOBuff;
  Tb *dCIOBuff, *eCIOBuff;
  Ta *zCIOBuff;
  char jobz = 'N';
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::steqr( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevd_external( FLA_EVD_WITHOUT_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<Tb>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<Tb>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<Ta>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "steqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "steqr(): Success\n");
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );

}

template< typename T >
void stedc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  T *dCPPIOBuff, *eCPPIOBuff, *zCPPIOBuff;
  T *dCIOBuff, *eCIOBuff, *zCIOBuff;
  char jobz = 'V';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::stedc( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevdd_external( FLA_EVD_WITH_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stedc(): Failure Diff = %E\n", diff);
  }else{
    printf( "stedc(): Success\n");
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );
}

template< typename Ta, typename Tb >
void stedc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  Tb *dCPPIOBuff, *eCPPIOBuff;
  Ta *zCPPIOBuff;
  Tb *dCIOBuff, *eCIOBuff;
  Ta *zCIOBuff;
  char jobz = 'V';
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::stedc( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevdd_external( FLA_EVD_WITH_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<Tb>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<Tb>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<Ta>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stedc(): Failure Diff = %E\n", diff);
  }else{
    printf( "stedc(): Success\n");
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );
}

//stedc_test();
//stemr_test();

template< typename T >
void syev_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, lCIOObj;
  T *aCPPIOBuff, *wCPPOBuff;
  T *aCIOBuff, *wCOBuff, *lCIOBuff;
  char jobz = 'N';
  char uplo = 'U';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  aCIOBuff  = new T[n*n];
  lCIOBuff  = new T[n*n];
  for(int j=0; j<n*n; j++)
  {
	  aCIOBuff[j] = 0;
	  lCIOBuff[j] = 0;
}

  //Call CPP function
  libflame::syev( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &lCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( lCIOBuff, 1, n, &lCIOObj );

  //Call C function
 // FLA_Error FLA_Hevdr_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l, FLA_Obj Z )
  FLA_Hevdr_external( FLA_EVD_WITHOUT_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj, lCIOObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  //diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syev(): Failure Diff = %E\n", diff);
  }else{
    printf( "syev(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &lCIOObj );
}

template< typename Ta, typename Tb >
void heev_test()
{
	#if 0
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, lCIOObj;
  T *aCPPIOBuff, *wCPPOBuff;
  T *aCIOBuff, *wCOBuff, *lCIOBuff;
  char jobz = 'V';
  char uplo = 'U';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  aCIOBuff  = new T[n*n];

  //Call CPP function
  libflame::syev( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &lCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( lCIOBuff, 1, n, &lCIOObj );

  //Call C function
 // FLA_Error FLA_Hevdr_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l, FLA_Obj Z )
  FLA_Hevdr_external( FLA_EVD_WITH_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj, lCIOObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syev(): Failure Diff = %E\n", diff);
  }else{
    printf( "syev(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &lCIOObj );
  #endif
}

template< typename T >
void syevd_test()
{
  int n = 64;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, lCIOObj;
  T *aCPPIOBuff, *wCPPOBuff;
  T *aCIOBuff, *wCOBuff, *lCIOBuff;
  char jobz = 'N';
  char uplo = 'U';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  aCIOBuff  = new T[n*n];
  lCIOBuff  = new T[n*n];

  //Call CPP function
  libflame::syevd( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &lCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( lCIOBuff, 1, n, &lCIOObj );

  //Call C function
  FLA_Hevdd_external( FLA_EVD_WITHOUT_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  //diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );
FILE *fp = fopen("test/in.txt", "a+");
print(fp, n, wCOBuff, wCPPOBuff);
fclose(fp);

  if(diff != 0.0)
  {
    printf( "syevd(): Failure Diff = %E\n", diff);
  }else{
    printf( "syevd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &lCIOObj );
}

void potrf_testall_variants(){
  potrf_test<float>();
  potrf_test<double>();
  potrf_test<lapack_complex_float>();
  potrf_test<lapack_complex_double>();
}

void potf2_testall_variants(){
  potf2_test<float>();
  potf2_test<double>();
  potf2_test<lapack_complex_float>();
  potf2_test<lapack_complex_double>();
}
void getrf_testall_variants(){
  getrf_test<float>();
  getrf_test<double>();
  getrf_test<lapack_complex_float>();
  getrf_test<lapack_complex_double>();
}

void getf2_testall_variants(){
  getf2_test<float>();
  getf2_test<double>();
  getf2_test<lapack_complex_float>();
  getf2_test<lapack_complex_double>();
}

void geqrf_testall_variants(){
  geqrf_test<float>();
  geqrf_test<double>();
  geqrf_test<lapack_complex_float>();
  geqrf_test<lapack_complex_double>();
}

void geqr2_testall_variants(){
  geqr2_test<float>();
  geqr2_test<double>();
  geqr2_test<lapack_complex_float>();
  geqr2_test<lapack_complex_double>();
}

void geqpf_testall_variants(){
  geqpf_test<float>();
  geqpf_test<double>();
  geqpf_test<lapack_complex_float, float>();
  geqpf_test<lapack_complex_double, double>();
}

void geqp3_testall_variants(){
  geqp3_test<float>();
  geqp3_test<double>();
  geqp3_test<lapack_complex_float, float>();
  geqp3_test<lapack_complex_double, double>();
}

void gelqf_testall_variants(){
  gelqf_test<float>();
  gelqf_test<double>();
  gelqf_test<lapack_complex_float>();
  gelqf_test<lapack_complex_double>();
}

void gelq2_testall_variants(){
  gelq2_test<float>();
  gelq2_test<double>();
  gelq2_test<lapack_complex_float>();
  gelq2_test<lapack_complex_double>();
}

void gelsd_testall_variants(){
  gelsd_test<float>();
  gelsd_test<double>();
  gelsd_test<lapack_complex_float, float>();
  gelsd_test<lapack_complex_double, double>();
}
void gelss_testall_variants(){
  gelss_test<float>();
  gelss_test<double>();
  gelss_test<lapack_complex_float, float>();
  gelss_test<lapack_complex_double, double>();
}

void lauum_testall_variants(){
  lauum_test<float>();
  lauum_test<double>();
  lauum_test<lapack_complex_float>();
  lauum_test<lapack_complex_double>();
}

void lauu2_testall_variants(){
  lauu2_test<float>();
  lauu2_test<double>();
  lauu2_test<lapack_complex_float>();
  lauu2_test<lapack_complex_double>();
}

void potri_testall_variants(){
  potri_test<float>();
  potri_test<double>();
  potri_test<lapack_complex_float>();
  potri_test<lapack_complex_double>();
}

void trtri_testall_variants(){
  trtri_test<float>();
  trtri_test<double>();
  trtri_test<lapack_complex_float>();
  trtri_test<lapack_complex_double>();
}

void trti2_testall_variants(){
  trti2_test<float>();
  trti2_test<double>();
  trti2_test<lapack_complex_float>();
  trti2_test<lapack_complex_double>();
}

void trsyl_testall_variants(){
  trsyl_test<float>();
  trsyl_test<double>();
  trsyl_test<lapack_complex_float, float>();
  trsyl_test<lapack_complex_double, double>();
}

void gehrd_testall_variants(){
  gehrd_test<float>();
  gehrd_test<double>();
  gehrd_test<lapack_complex_float>();
  gehrd_test<lapack_complex_double>();
}

void gehd2_testall_variants(){
  gehd2_test<float>();
  gehd2_test<double>();
  gehd2_test<lapack_complex_float>();
  gehd2_test<lapack_complex_double>();
}

void sytrd_testall_variants(){
  sytrd_test<float>();
  sytrd_test<double>();
}

void hetrd_testall_variants(){
  hetrd_test<lapack_complex_float, float>();
  hetrd_test<lapack_complex_double, double>();
}

void sytd2_testall_variants(){
  sytd2_test<float>();
  sytd2_test<double>();
}

void hetd2_testall_variants(){
  hetd2_test<lapack_complex_float, float>();
  hetd2_test<lapack_complex_double, double>();
}

void gebrd_testall_variants(){
  gebrd_test<float>();
  gebrd_test<double>();
  gebrd_test<lapack_complex_float, float>();
  gebrd_test<lapack_complex_double, double>();
}

void gebd2_testall_variants(){
  gebd2_test<float>();
  gebd2_test<double>();
  gebd2_test<lapack_complex_float, float>();
  gebd2_test<lapack_complex_double, double>();
}

void sygst_testall_variants(){
  sygst_test<float>();
  sygst_test<double>();
}

void hegst_testall_variants(){
  hegst_test<lapack_complex_float>();
  hegst_test<lapack_complex_double>();
}

void sygs2_testall_variants(){
  sygs2_test<float>();
  sygs2_test<double>();
}

void hegs2_testall_variants(){
  hegs2_test<lapack_complex_float>();
  hegs2_test<lapack_complex_double>();
}

void larft_testall_variants(){
  larft_test<float>();
  larft_test<double>();
  larft_test<lapack_complex_float>();
  larft_test<lapack_complex_double>();
}

void larfg_testall_variants(){
  larfg_test<float>();
  larfg_test<double>();
  larfg_test<lapack_complex_float>();
  larfg_test<lapack_complex_double>();
}

void larfgp_testall_variants(){
  larfgp_test<float>();
  larfgp_test<double>();
  larfgp_test<lapack_complex_float>();
  larfgp_test<lapack_complex_double>();
}

void orgqr_testall_variants(){
  orgqr_test<float>();
  orgqr_test<double>();
}
void ungqr_testall_variants(){
  ungqr_test<lapack_complex_float>();
  ungqr_test<lapack_complex_double>();
}
void ormqr_testall_variants(){
  ormqr_test<float>();
  ormqr_test<double>();
}
void unmqr_testall_variants(){
  unmqr_test<lapack_complex_float>();
  unmqr_test<lapack_complex_double>();
}
void orm2r_testall_variants(){
  orm2r_test<float>();
  orm2r_test<double>();
}
void unm2r_testall_variants(){
  unm2r_test<lapack_complex_float>();
  unm2r_test<lapack_complex_double>();
}

void orglq_testall_variants(){
  orglq_test<float>();
  orglq_test<double>();
  orglq_test<lapack_complex_float>();
  orglq_test<lapack_complex_double>();
}
void ormlq_testall_variants(){
  ormlq_test<float>();
  ormlq_test<double>();
}
void unmlq_testall_variants(){
  unmlq_test<lapack_complex_float>();
  unmlq_test<lapack_complex_double>();
}
  //orml2
  //unml2
void orgtr_testall_variants(){
  orgtr_test<float>();
  orgtr_test<double>();
}
void ungtr_testall_variants(){
  ungtr_test<lapack_complex_float, float>();
  ungtr_test<lapack_complex_double, double >();
}
void ormtr_testall_variants(){
  ormtr_test<float>();
  ormtr_test<double>();
}
void unmtr_testall_variants(){
  unmtr_test<lapack_complex_float, float>();
  unmtr_test<lapack_complex_double, double >();
}
void orgbr_testall_variants(){
  orgbr_test<float>();
  orgbr_test<double>();
}
void ungbr_testall_variants(){
  ungbr_test<lapack_complex_float, float>();
  ungbr_test<lapack_complex_double, double >();
}
void ormbr_testall_variants(){
  ormbr_test<float>();
  ormbr_test<double>();
}
void unmbr_testall_variants(){
  unmbr_test<lapack_complex_float, float>();
  unmbr_test<lapack_complex_double, double >();
}
void steqr_testall_variants(){
  steqr_test<float>();
  steqr_test<double>();
  steqr_test<lapack_complex_float, float>();
  steqr_test<lapack_complex_double, double >();
}
void stedc_testall_variants(){
  stedc_test<float>();
  stedc_test<double>();
  stedc_test<lapack_complex_float, float>();
  stedc_test<lapack_complex_double, double >();
}
//stedc_testall_variants();
//stemr_testall_variants();
void syev_testall_variants(){
  syev_test<float>();
  syev_test<double>();
}
void heev_testall_variants(){
  heev_test<lapack_complex_float, float>();
  heev_test<lapack_complex_double, double >();
}
void syevd_testall_variants(){
  syevd_test<float>();
  syevd_test<double>();
}
//heevd_testall_variants();
//syevr_testall_variants();


//void unglq_testall_variants(){
//  unglq_test<float>();
//  unglq_test<double>();
//  unglq_test<lapack_complex_float>();
//  unglq_test<lapack_complex_double>();
//}

//#define Test(fnName)\
// test_ ## fnName ## (float)();
//
 //test_ ## fnName <double>();\
 //test_ ## fnName <lapack_complex_float>();\
 //test_ ## fnName <lapack_complex_double>();

int main(int argc, char *argv[])
{
  // steqr_testall_variants();//pass
  // stedc_testall_variants();//pass
  //stemr_testall_variants();
  //syev_testall_variants(); //fail
  //heev_testall_variants();
  //syevd_testall_variants(); // cout is zero
  //heevd_testall_variants();
  //syevr_testall_variants();
  //heevr
  //heevr
  //bdsqr
  //bdsdc
  //gesvd
  //gesdd
  //laswp
  //laset



  //larft_testall_variants(); //larft not found
  //larfgp_testall_variants(); //slarfg_ not defined
  //orm2r_testall_variants();//sorm2r_ not defined
  //unm2r_testall_variants();//sunm2r_ not defined
  //orglq_testall_variants(); //Pending
  //unglq_testall_variants(); //Pending
  //ormlq_testall_variants(); //pass
  //unmlq_testall_variants(); //pass
  //orml2_testall_variants(); //pass
  //unml2_testall_variants(); //pass
  //trsyl_testall_variants(); //pass
  //gehd2_testall_variants(); //pass
#if 1
  potrf_testall_variants(); //pass
  potf2_testall_variants(); //pass
  getrf_testall_variants(); //pivot mismatch
  getf2_testall_variants();//pivot mismatch
  geqrf_testall_variants(); //pass
  geqr2_testall_variants(); //pass
  geqpf_testall_variants(); //pass  //m,n>128 fails
  geqp3_testall_variants(); //pass  //m,n>128 fails
  gelqf_testall_variants(); //pass
  gelq2_testall_variants(); //pass
  gelsd_testall_variants(); //pass
  gelss_testall_variants(); //pass
  lauum_testall_variants(); //pass
  lauu2_testall_variants(); //pass
  potri_testall_variants(); //pass
  trtri_testall_variants(); //pass //m>128 fails
  trti2_testall_variants(); //pass //m>128 fails
  trsyl_testall_variants(); //pass
  gehrd_testall_variants();//pass
  gehd2_testall_variants(); //pass //tau buff  memset to 0
  sytrd_testall_variants(); //pass
  hetrd_testall_variants();//pass
  sytd2_testall_variants(); //pass
  //hetd2_testall_variants(); //fails //in/out buff failure--after 3 decimal point
  gebrd_testall_variants(); //pass
  //gebd2_testall_variants(); //fails //alll buff mismatch-- after few decimal point
  //    /*Testing to be done*/
  sygst_testall_variants();//pass
  hegst_testall_variants();//pass //m>128 fails for complex float
  //sygs2_testall_variants();//pass //m>64 fails for complex float
  hegs2_testall_variants();//pass //m>64 fails for complex float
  larft_testall_variants(); //larft not found
  larfg_testall_variants(); //pass
  //larfgp_testall_variants(); //slarfg_ not defined
  orgqr_testall_variants(); //pass
  //ungqr_testall_variants(); //fails
  ormqr_testall_variants(); //pass
  unmqr_testall_variants(); //pass
  //orm2r
  //unm2r
  //orglq_testall_variants();//implemenation pending
  //unglq_testall_variants();//implemenation pending
  ormlq_testall_variants(); //pass
  unmlq_testall_variants(); //pass
  //orml2
  //unml2
  orgtr_testall_variants(); //pass
  ungtr_testall_variants(); //pass
  ormtr_testall_variants(); //pass
  unmtr_testall_variants(); //pass
  orgbr_testall_variants(); //pass
  ungbr_testall_variants(); //pass
  ormbr_testall_variants(); //pass
  unmbr_testall_variants(); //pass
  //steqr_testall_variants(); //pass
  //stedc_testall_variants(); //pass
  //stemr
  //syev
  //heev
  //syevd
  //heevd
  //syevr
  //heevr
  //heevr
  //bdsqr
  //bdsdc
  //gesvd
  //gesdd
  //laswp
  //laset
  #endif
}


