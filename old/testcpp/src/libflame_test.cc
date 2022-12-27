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

FLA_Error larft_C(FLA_Direct direct, FLA_Direct storev, FLA_Obj A, FLA_Obj tauObj, FLA_Obj T)
{
  int          info = 0;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  datatype = FLA_Obj_datatype( A );

  int m_A      = FLA_Obj_length( A );
  int n_A      = FLA_Obj_width( A );
  int cs_A     = FLA_Obj_col_stride( A );
  int ldt      = FLA_Obj_length( T );
  char char_direct,  char_storev;
  FLA_Param_map_flame_to_netlib_direct( direct, &char_direct );
  FLA_Param_map_flame_to_netlib_storev( storev, &char_storev );
  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_T    = ( float * ) FLA_FLOAT_PTR( T );
      float *buff_tau    = ( float * ) FLA_FLOAT_PTR( tauObj );
      slarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_T    = ( double * ) FLA_DOUBLE_PTR( T );
      double *buff_tau    = ( double * ) FLA_DOUBLE_PTR( tauObj );
      dlarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_T    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( T );
      lapack_complex_float *buff_tau    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( tauObj );
      clarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_T    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( T );
      lapack_complex_double *buff_tau    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( tauObj );
      zlarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error laset_C( FLA_Uplo uplo, FLA_Obj alphaObj, FLA_Obj betaObj, FLA_Obj A)
{
  int          info = 0;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  datatype = FLA_Obj_datatype( A );
  int m_A      = FLA_Obj_length( A );
  int n_A      = FLA_Obj_width( A );
  int cs_A     = FLA_Obj_col_stride( A );
  char blas_uplo;
  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *alpha    = ( float * ) FLA_FLOAT_PTR( alphaObj );
      float *beta    = ( float * ) FLA_FLOAT_PTR( betaObj );
      slaset_( &blas_uplo, &m_A, &n_A, alpha, beta, buff_A, &cs_A );
      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *alpha    = ( double * ) FLA_DOUBLE_PTR( alphaObj );
      double *beta    = ( double * ) FLA_DOUBLE_PTR( betaObj );
      dlaset_( &blas_uplo, &m_A, &n_A, alpha, beta, buff_A, &cs_A );
      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *alpha    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( alphaObj );
      lapack_complex_float *beta    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( betaObj );
       claset_( &blas_uplo, &m_A, &n_A, alpha, beta, buff_A, &cs_A );
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *alpha    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( alphaObj );
      lapack_complex_double *beta    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( betaObj );
      zlaset_( &blas_uplo, &m_A, &n_A, alpha, beta, buff_A, &cs_A );
      break;
    }
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}


FLA_Error orglq_C( FLA_Obj A, FLA_Obj t, int k )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          cs_A;
  FLA_Obj work_obj;
  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );
  int lwork    = m_A * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );
  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
     sorglq_( &m_A, &n_A, &k, buff_A, &cs_A, buff_t, buff_work, &lwork, &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );
    dorglq_( &m_A, &n_A, &k, buff_A, &cs_A, buff_t, buff_work, &lwork, &info );

    break;
  }

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
    lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
    lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
    cunglq_( &m_A, &n_A, &k, buff_A, &cs_A, buff_t, buff_work, &lwork, &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
    lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
    lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
    zunglq_( &m_A, &n_A, &k, buff_A, &cs_A, buff_t, buff_work, &lwork, &info );

    break;
  }

  FLA_Obj_free( &work_obj );
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error unml2_C( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
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

       sormlq_( &blas_side,
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

     dormlq_( &blas_side,
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
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
      lapack_complex_float *buff_B    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( B );
      lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );

      cunmlq_( &blas_side,
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
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
      lapack_complex_double *buff_B    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( B );
      lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

      zunmlq_( &blas_side,
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
    printf( "potrf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "potrf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "potf2():%s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "potf2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );

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
    printf( "%s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else if(diffInt !=0){
    printf( "getrf(): Failure Pivot Buffer mismatach Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );

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
    printf( "%s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else if(diffInt !=0){
    printf( "getf2(): Failure:Pivot Buffer Mismatch Diff = %d\n", diffInt);
  }else{
    printf( "getf2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "geqrf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "geqrf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;

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
    printf( "%s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "geqr2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqr2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 1024;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  T *aCPPIOBuff, *aCIOBuff, *rworkRefBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = fla_min( m, n ) ;

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
    printf( "geqpf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = fla_min( m, n ) ;

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
    printf( "geqpf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 2048;
  int n = 1024;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  T *aCPPIOBuff, *aCIOBuff, *rworkRefBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = fla_min( m, n ) ;

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
    printf( "geqp3(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 128;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = fla_min( m, n ) ;

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
    printf( "geqp3(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff1111 = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;

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
    printf( "%s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "gelqf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelqf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;
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
    printf( "%s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else if(diffInt !=0){
    printf( "gelq2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelq2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;
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
    printf( "gelsd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else{
    printf( "gelsd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;
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
    printf( "gelsd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else{
    printf( "gelsd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;
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
    printf( "gelss(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else{
    printf( "gelss(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n ) ;
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
    printf( "gelss(): %s TEST FAIL\n" , __PRETTY_FUNCTION__) ;
  }else{
    printf( "gelss(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "lauum(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "lauum(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "lauu2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "lauu2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "potri(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "potri(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "trtri(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "trtri(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "trti2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "trti2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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

  int datatype = getDatatype<T>();

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
  double diff =  computeError<T>( m, m, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
   printf( "trsyl(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "trsyl(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
   printf( "trsyl(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "gehrd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gehrd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int iho = 1;
  int ilo = 1;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);

  tauCPPOBuff =  new T [n-1];
  tauCOBuff =  new T [n-1];
  for(int i=0; i < n-1; i++)
  {
    tauCPPOBuff[i] = 0;
    tauCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gehd2( LAPACK_COL_MAJOR, &n, &ilo, &iho, aCPPIOBuff, &n, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Hess_unb_external( aCIOObj, tauCOObj, 1, 1 );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gehd2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gehd2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "sytrd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "sytrd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "hetrd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "hetrd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "sytd2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "sytd2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
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

  //Call CPP function
  libflame::hetd2( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "hetd2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "hetd2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
void gebrd_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = fla_min( m, n );
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
    printf( "gebrd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gebrd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "gebrd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gebrd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<T>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebd2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gebd2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff);

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
    printf( "gebd2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gebd2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  char uplo  ='L';
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
    printf( "sygst(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "sygst(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "hegst(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "hegst(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int n = 2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 2 ; //1 or 2
  int itype_c = FLA_NO_INVERSE;

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
    printf( "sygs2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "sygst(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 2 ; //1 or 2
  int itype_c = FLA_NO_INVERSE;

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
    printf( "hegs2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "hegs2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int n = 64;//256;
  int k = 64;//128;
 srand (time(NULL));

 FLA_Init( );
 FLA_Obj vCIObj, tauCIObj, tCOObj;
 T *vCPPIBuff, *vCIBuff, *tCPPOBuff, *tCOBuff ;
 T *tauCPPIBuff, *tauCIBuff ;
 int datatype = getDatatype<T>();
 char direct  = 'F';
 char storev = 'C';
 int ldv = n;
 int ldt = k;
 //Allocate and initialize buffers for C and CPP functions with random values
 allocate_init_buffer(vCPPIBuff, vCIBuff, ldv*k);
 tauCPPIBuff =  new T [k];
 tauCIBuff =  new T [k];
 tCPPOBuff =  new T [ldt*k];
 tCOBuff =  new T [ldt*k];
 for(int i =0; i <ldt*k; i++)
 {
    tCPPOBuff[i] = 0;
    tCOBuff[i] = 0;
 }

 //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, vCPPIBuff, &n, tauCPPIBuff );
  libflame::larft( LAPACK_COL_MAJOR, &direct, &storev, &n, &k, vCPPIBuff, &ldv, tauCPPIBuff, tCPPOBuff, &ldt );

 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, ldv, k, &vCIObj );
 FLA_Obj_create_without_buffer( datatype, k, 1, &tauCIObj );
 FLA_Obj_create_without_buffer( datatype, ldt, k, &tCOObj );
 FLA_Obj_attach_buffer( vCIBuff, 1, ldv, &vCIObj );
 FLA_Obj_attach_buffer( tauCIBuff, 1, k, &tauCIObj );
 FLA_Obj_attach_buffer( tCOBuff, 1, ldt, &tCOObj );

  //Call C function
  FLA_QR_blk_external( vCIObj, tauCIObj );
  larft_C( FLA_FORWARD, FLA_COLUMNWISE, vCIObj, tauCIObj, tCOObj );
  double diff =  computeError<T>( ldt, k, tCOBuff, tCPPOBuff );

  if(diff != 0.0)
  {
    printf( "larft(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "larft(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete vCPPIBuff ;
  delete tauCPPIBuff ;
  delete tCPPOBuff;
  FLA_Obj_free( &vCIObj );
  FLA_Obj_free( &tauCIObj );
  FLA_Obj_free( &tCOObj );
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
    printf( "larfg(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }
  else{
    printf( "larfg(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete xCPPIOBuff ;
  FLA_Obj_free( &xCIOObj );
}

template< typename T >
void larfgp_test()
{
  int n = 1024;
  srand (time(NULL));

  FLA_Init( );
  T alphaC = 1, alphaCPP =1;
  T tauC = 0;
  T tauCPP = 0;
  FLA_Obj xCIOObj, alphaCIOObj, tauCOObj;
  T *xCPPIOBuff, *xCIOBuff;
  int datatype = getDatatype<T>();
  int incxCPP = 5;
  int incx = 5;
  int size  = (1+(n-2)*abs(incx));
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(xCPPIOBuff, xCIOBuff, size);

  //Call CPP function
  libflame::larfg( &n, &alphaCPP, xCPPIOBuff, &incxCPP, &tauCPP);

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
    printf( "larfgp(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }
  else{
    printf( "larfgp(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
 int min_m_n = fla_min( m, n );
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
    printf( "orgqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orgqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int k = 256;
  int min_m_n = fla_min( m, n );
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
    printf( "ungqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ungqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "ormqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "unmqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unmqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 512;
  int n = 256;
  int k = 128;
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
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orm2r(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orm2r(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 64;
  int n = 64;
  int k = 64;
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
  FLA_Obj_attach_buffer( cCIOBuff, 1, n, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unm2r(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unm2r(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int m = 256;
  int n = 512;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  int min_m_n = fla_min(m,n);
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
  libflame::orglq( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  orglq_C( aCIOObj, tauCOObj, k );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orglq(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orglq(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void unglq_test()
{
  int m = 256;
  int n = 512;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  int min_m_n = fla_min(m,n);

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
  libflame::unglq( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  orglq_C( aCIOObj, tauCOObj, k );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "unglq(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unglq(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "ormlq(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormlq(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "unmlq(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unmlq(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void orml2_test()
{
  int m = 512;
  int n = 256;
  int k = 128;
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
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];
  cCIOBuff =  new T [m*n];
  cCPPIOBuff =  new T [m*n];
  for (int i = 0; i < m*n; i++)
  {
    cCPPIOBuff[i] = 0;
    cCIOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &k, &n, aCPPIOBuff, &k, tauCPPOBuff );
  libflame::orml2( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &k, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, k, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, k, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  unml2_C( FLA_RIGHT, FLA_TRANSPOSE, FLA_ROWWISE, aCIOObj, tauCOObj, cCIOObj );

  double diff =  computeError<T>( k, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orml2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orml2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );

}

template< typename T >
void unml2_test()
{
  int m = 256;
  int n = 512;
  int k = 128;
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
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  cCIOBuff =  new T [m*n];
  cCPPIOBuff =  new T [m*n];
  for (int i = 0; i < m*n; i++)
  {
    cCPPIOBuff[i] = 0;
    cCIOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &k, &n, aCPPIOBuff, &k, tauCPPOBuff );
  libflame::unml2( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &k, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, k, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, k, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  unml2_C( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_ROWWISE, aCIOObj, tauCOObj, cCIOObj );

  double diff =  computeError<T>( k, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unml2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unml2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "orgtr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orgtr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "ungtr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ungtr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "ormtr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormtr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "unmtr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unmtr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "orgbr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orgbr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "ungbr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ungbr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "ormbr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormbr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
  int min_m_n = fla_min( m, n );
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
    printf( "unmbr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "unmbr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "steqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "steqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "steqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "steqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "stedc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stedc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
    printf( "stedc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stedc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
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
void stemr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj, eCIOObj;
  T *aCPPIOBuff, *wCPPOBuff, *zCPPOBuff, *eCPPIOBuff;
  int *isuppzCPPOBuff;
  T *aCIOBuff, *wCOBuff, *zCOBuff, *eCIOBuff;
  char jobz = 'N'; //N, V
  char range = 'A'; //A, V, I
  int datatype = getDatatype<T>();
  int v1 = 0, vu = 0;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  int tryrac = 1;
  int nzc = n;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  zCOBuff  = new T[n*m];
  zCPPOBuff = new T[n*m];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*m ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::stemr( LAPACK_COL_MAJOR, &jobz, &range, &n, aCPPIOBuff, eCPPIOBuff, &v1, &vu, &il, &iu,
                   &m, wCPPOBuff, zCPPOBuff, &ldz, &nzc, isuppzCPPOBuff, &tryrac);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, m, &zCOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &eCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, n, &eCIOObj );

  //Call C function
  FLA_Tevdr_external( FLA_EVD_WITHOUT_VECTORS, aCIOObj, eCIOObj, wCOObj, zCOObj );

  double diff =  computeError<T>( n, 1, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<T>( n, m, zCOBuff, zCPPOBuff );
  diff +=  computeError<T>( n, 1, eCIOBuff, eCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stemr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stemr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete eCPPIOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
  FLA_Obj_free( &eCIOObj );
}

template< typename Ta, typename Tb >
void stemr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj, eCIOObj;
   Ta *zCPPOBuff ;
  Tb *wCPPOBuff, *aCPPIOBuff, *eCPPIOBuff;
  int *isuppzCPPOBuff;
  Ta *zCOBuff;
  Tb *wCOBuff, *aCIOBuff, *eCIOBuff;
  char jobz = 'N'; //N, V
  char range = 'A'; //A, V, I
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  int v1 = 0, vu = 0;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  int tryrac = 1;
  int nzc = n;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];
  zCOBuff  = new Ta[n*m];
  zCPPOBuff = new Ta[n*m];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*m ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::stemr( LAPACK_COL_MAJOR, &jobz, &range, &n, aCPPIOBuff, eCPPIOBuff, &v1, &vu, &il, &iu,
                   &m, wCPPOBuff, zCPPOBuff, &ldz, &nzc, isuppzCPPOBuff, &tryrac);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, m, &zCOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &eCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, n, &eCIOObj );

  //Call C function
  FLA_Tevdr_external( FLA_EVD_WITHOUT_VECTORS, aCIOObj, eCIOObj, wCOObj, zCOObj );

  double diff =  computeError<Tb>( n, 1, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<Ta>( n, m, zCOBuff, zCPPOBuff );
  diff +=  computeError<Tb>( n, 1, eCIOBuff, eCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stemr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stemr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete eCPPIOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
  FLA_Obj_free( &eCIOObj );
}



template< typename T >
void syev_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj;
  T *aCPPIOBuff, *wCPPOBuff;
  T *aCIOBuff, *wCOBuff;
  char jobz = 'N';
  char uplo = 'U';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];

  //Call CPP function
  libflame::syev( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );

  //Call C function
  FLA_Hevd_external( FLA_EVD_WITHOUT_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syev(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "syev(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
}

template< typename Ta, typename Tb >
void heev_test()
{
  int n = 2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj;
  Ta *aCPPIOBuff, *aCIOBuff;
  Tb *wCPPOBuff, *wCOBuff;
  char jobz = 'N';
  char uplo = 'U';
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];

  //Call CPP function
  libflame::heev( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );

  //Call C function
  FLA_Hevd_external( FLA_EVD_WITHOUT_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj );

  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "heev(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "heev(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
}

template< typename T >
void syevd_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj;
  T *aCPPIOBuff, *wCPPOBuff;
  T *aCIOBuff, *wCOBuff;
  char jobz = 'N';
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];

  //Call CPP function
  libflame::syevd( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );

  //Call C function
  FLA_Hevdd_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, aCIOObj, wCOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
 // diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syevd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "syevd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
}

template< typename Ta, typename Tb >
void heevd_test()
{
  int n = 2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj;
  Ta *aCPPIOBuff, *aCIOBuff;
  Tb *wCPPOBuff, *wCOBuff;
  char jobz = 'N';
  char uplo = 'U';
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];

  //Call CPP function
  libflame::heevd( LAPACK_COL_MAJOR, &jobz, &uplo, &n, aCPPIOBuff, &n, wCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );

  //Call C function
  FLA_Hevdd_external( FLA_EVD_WITHOUT_VECTORS, FLA_UPPER_TRIANGULAR, aCIOObj, wCOObj );

  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );

  if(diff != 0.0)
  {
    printf( "heevd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "heevd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
}

template< typename T >
void syevr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj;
  T *aCPPIOBuff, *wCPPOBuff, *zCPPOBuff;
  int *isuppzCPPOBuff;
  T *aCIOBuff, *wCOBuff, *zCOBuff;
  char jobz = 'N'; //N, V
  char uplo = 'l';
  char range = 'A'; //A, V, I
  int datatype = getDatatype<T>();
  T v1 = 0, vu = 0, abstol = FLA_MACH_SFMIN;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  zCOBuff  = new T[n*n];
  zCPPOBuff = new T[n*n];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*n ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::syevr( LAPACK_COL_MAJOR, &jobz, &range, &uplo, &n, aCPPIOBuff, &n, &v1, &vu, &il, &iu,
                   &abstol, &m, wCPPOBuff, zCPPOBuff, &ldz, isuppzCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );

  //Call C function
  FLA_Hevdr_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, aCIOObj, wCOObj, zCOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<T>( n, n, zCOBuff, zCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syevr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "syevr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
}

template< typename Ta, typename Tb >
void heevr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj;
  Ta *aCPPIOBuff, *zCPPOBuff;
  Tb  *wCPPOBuff;
  int *isuppzCPPOBuff;
  Ta *aCIOBuff, *zCOBuff;
  Tb *wCOBuff;
  char jobz = 'N'; //N, V
  char uplo = 'l';
  char range = 'A'; //A, V, I
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  Tb v1 = 0, vu = 0, abstol = FLA_MACH_SFMIN;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];
  zCOBuff  = new Ta[n*n];
  zCPPOBuff = new Ta[n*n];
  isuppzCPPOBuff = new int [isuppzDim];

  for(int i =0; i < n*n ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::heevr( LAPACK_COL_MAJOR, &jobz, &range, &uplo, &n, aCPPIOBuff, &n, &v1, &vu, &il, &iu,
                   &abstol, &m, wCPPOBuff, zCPPOBuff, &ldz, isuppzCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &zCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );

  //Call C function
  FLA_Hevdr_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, aCIOObj, wCOObj, zCOObj );

  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<Ta>( n, n, zCOBuff, zCPPOBuff );

  if(diff != 0.0)
  {
    printf( "heevr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "heevr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
}

template< typename T >
void bdsqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  T *dCPPIOBuff, *dCIOBuff;
  T *eCPPIOBuff, *eCIOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPIBuff, *vtCIBuff;
  T *cCPPOBuff;
  int datatype = getDatatype<T>();
  char uplo = 'L';
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new T[n*n];
  uCOBuff = new T[n*n];
  cCPPOBuff = new T[n*n];
  int ncvt = n, nru = n, ncc = n, ldc = n;
  for(int i =0; i < n*n ; i++)
  {
    uCOBuff[i] = 0;
    uCPPOBuff[i] = 0;
  }


  //Call CPP function
  libflame::bdsqr( LAPACK_COL_MAJOR, &uplo, &n, &ncvt, &nru, &ncc, dCPPIOBuff, eCPPIOBuff, vtCPPIBuff, &n,  uCPPOBuff, &n, cCPPOBuff, &ldc);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCIObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, n, &uCOObj );
  FLA_Obj_attach_buffer( vtCIBuff, 1, n, &vtCIObj );

  //Call C function
  FLA_Bsvd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete uCPPOBuff ;
  delete vtCPPIBuff ;
  delete cCPPOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCIObj );
}

template< typename Ta, typename Tb >
void bdsqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  Tb *dCPPIOBuff, *dCIOBuff;
  Tb *eCPPIOBuff, *eCIOBuff;
  Ta *uCPPOBuff, *uCOBuff;
  Ta *vtCPPIBuff, *vtCIBuff;
  Ta *cCPPOBuff;
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  char uplo = 'L';
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new Ta[n*n];
  uCOBuff = new Ta[n*n];
  cCPPOBuff = new Ta[n*n];
  int ncvt = n, nru = n, ncc = n, ldc = n;
  //Call CPP function
  libflame::bdsqr( LAPACK_COL_MAJOR, &uplo, &n, &ncvt, &nru, &ncc, dCPPIOBuff, eCPPIOBuff, vtCPPIBuff, &n,  uCPPOBuff, &n, cCPPOBuff, &ldc);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &uCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &vtCIObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, n, &uCOObj );
  FLA_Obj_attach_buffer( vtCIBuff, 1, n, &vtCIObj );

  //Call C function
  FLA_Bsvd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<Tb>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<Tb>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<Ta>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<Ta>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete uCPPOBuff ;
  delete vtCPPIBuff ;
  delete cCPPOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCIObj );
}

template< typename T >
void bdsdc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  T *dCPPIOBuff, *dCIOBuff;
  T *eCPPIOBuff, *eCIOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPIBuff, *vtCIBuff;
  T *qCPPOBuff;
  int *iqCPPOBuff;
  int SMLSIZ =25;
  int LDQ = n*(11 + 2*SMLSIZ + 8*(log2(n/(SMLSIZ+1))));
  int  LDIQ = n*(3 + 3*(log2(n/(SMLSIZ+1))));
  int datatype = getDatatype<T>();
  char uplo = 'L';
  char compq = 'I'; //N, P, I
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new T[n*n];
  uCOBuff = new T[n*n];
  qCPPOBuff = new T[LDQ];
  iqCPPOBuff = new int[LDIQ];

  //Call CPP function
  libflame::bdsdc( LAPACK_COL_MAJOR, &uplo, &compq, &n, dCPPIOBuff, eCPPIOBuff, uCPPOBuff, &n, vtCPPIBuff, &n, qCPPOBuff, iqCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCIObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, n, &uCOObj );
  FLA_Obj_attach_buffer( vtCIBuff, 1, n, &vtCIObj );

  //Call C function
  FLA_Bsvdd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsdc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsdc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete uCPPOBuff ;
  delete vtCPPIBuff ;
  delete qCPPOBuff ;
  delete iqCPPOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCIObj );
}


template< typename T >
void gesvd_test()
{
  int m = 64;
  int n = 64;
  int min_m_n = fla_min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  T *aCPPIOBuff, *aCIOBuff;
  T *sCPPOBuff, *sCOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPOBuff, *vtCOBuff;
  T *superb;

  int datatype = getDatatype<T>();
  char jobu = 'A'; //A , S, O, N
  char jobv = 'A'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new T[min_m_n];
  sCOBuff = new T[min_m_n];
  uCPPOBuff = new T[m*m];
  uCOBuff = new T[m*m];
  vtCPPOBuff = new T[m*m];
  vtCOBuff = new T[m*m];
  superb = new T[min_m_n];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gesvd( LAPACK_COL_MAJOR, &jobu, &jobv, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n, superb );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatype, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svd_external( FLA_SVD_VECTORS_ALL, FLA_SVD_VECTORS_ALL, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<T>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gesvd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gesvd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}

template< typename Ta, typename Tb >
void gesvd_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = fla_min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  Ta *aCPPIOBuff, *aCIOBuff;
  Tb *sCPPOBuff, *sCOBuff;
  Ta *uCPPOBuff, *uCOBuff;
  Ta *vtCPPOBuff, *vtCOBuff;
  Tb *superb;

  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  char jobu = 'A'; //A , S, O, N
  char jobv = 'A'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new Tb[min_m_n];
  sCOBuff = new Tb[min_m_n];
  uCPPOBuff = new Ta[m*m];
  uCOBuff = new Ta[m*m];
  vtCPPOBuff = new Ta[m*m];
  vtCOBuff = new Ta[m*m];
  superb = new Tb[min_m_n];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gesvd( LAPACK_COL_MAJOR, &jobu, &jobv, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n, superb );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svd_external( FLA_SVD_VECTORS_ALL, FLA_SVD_VECTORS_ALL, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<Ta>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<Ta>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gesvd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gesvd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}

template< typename T >
void gesdd_test()
{
  int m = 64;
  int n = 64;
  int min_m_n = fla_min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  T *aCPPIOBuff, *aCIOBuff;
  T *sCPPOBuff, *sCOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPOBuff, *vtCOBuff;

  int datatype = getDatatype<T>();
  char jobz = 'O'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new T[min_m_n];
  sCOBuff = new T[min_m_n];
  uCPPOBuff = new T[m*m];
  uCOBuff = new T[m*m];
  vtCPPOBuff = new T[m*m];
  vtCOBuff = new T[m*m];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //int gesdd( int matrix_layout, char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt )
  //Call CPP function
  libflame::gesdd( LAPACK_COL_MAJOR, &jobz, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatype, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svdd_external( FLA_SVD_VECTORS_MIN_OVERWRITE, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<T>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gessd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gessd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}

template< typename Ta, typename Tb >
void gesdd_test()
{
  int m = 64;
  int n = 64;
  int min_m_n = fla_min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  Ta *aCPPIOBuff, *aCIOBuff;
  Tb *sCPPOBuff, *sCOBuff;
  Ta *uCPPOBuff, *uCOBuff;
  Ta *vtCPPOBuff, *vtCOBuff;

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  char jobz = 'O'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new Tb[min_m_n];
  sCOBuff = new Tb[min_m_n];
  uCPPOBuff = new Ta[m*m];
  uCOBuff = new Ta[m*m];
  vtCPPOBuff = new Ta[m*m];
  vtCOBuff = new Ta[m*m];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //int gesdd( int matrix_layout, char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt )
  //Call CPP function
  libflame::gesdd( LAPACK_COL_MAJOR, &jobz, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatype, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svdd_external( FLA_SVD_VECTORS_MIN_OVERWRITE, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<Ta>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<Ta>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gessd(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "gessd(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}


template< typename T >
void laswp_test()
{
  int n = 1024;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, ipivIObj;
  T *aCPPIOBuff, *aCIOBuff;
  int *ipivCPPIBuff, *ipivCIBuff;
  int k1 = 1;
  int k2 =  128;
  int incx = 1;
  int pivDim = (k1+(k2-k1)*abs(incx));
  int datatype = getDatatype<T>();
  int datatypeTb = getDatatype<int>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(ipivCPPIBuff, ipivCIBuff, pivDim);

  for (int i = 0; i < pivDim; i++ )
  {
    ipivCPPIBuff[ i ] = ipivCIBuff[ i ] + i + 1;
  }

  //Call CPP function
  libflame::laswp( LAPACK_COL_MAJOR, &n, aCPPIOBuff, &n, &k1, &pivDim, ipivCPPIBuff, &incx );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, pivDim, 1, &ipivIObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( ipivCIBuff, 1, pivDim, &ipivIObj );

  //Call C function
  FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, ipivIObj, aCIOObj  );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "laswp(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "laswp(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete ipivCPPIBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &ipivIObj );
}

template< typename T >
void laset_test()
{
  int m = 1024;
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, alphaObj, betaObj;
  T *aCPPIOBuff;
  T *aCIOBuff;
  T alpha =1;
  T beta =1;
  char uplo = 'U';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  //Call CPP function
  libflame::laset( LAPACK_COL_MAJOR, &uplo, &m, &n, &alpha, &beta, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &alphaObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &betaObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( &alpha, 1, 1, &alphaObj );
  FLA_Obj_attach_buffer( &beta, 1, 1, &betaObj );

  //Call C function
  laset_C( FLA_UPPER_TRIANGULAR, alphaObj, betaObj, aCIOObj  );

  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "laset(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "laset(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  FLA_Obj_free( &aCIOObj );
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
}
void unglq_testall_variants(){
  unglq_test<lapack_complex_float>();
  unglq_test<lapack_complex_double>();
}
void ormlq_testall_variants(){
  ormlq_test<float>();
  ormlq_test<double>();
}
void unmlq_testall_variants(){
  unmlq_test<lapack_complex_float>();
  unmlq_test<lapack_complex_double>();
}
void orml2_testall_variants(){
  orml2_test<float>();
  orml2_test<double>();
}
void unml2_testall_variants(){
  unml2_test<lapack_complex_float>();
  unml2_test<lapack_complex_double>();
}
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
void stemr_testall_variants(){
  stemr_test<float>();
  stemr_test<double>();
  stemr_test<lapack_complex_float, float>();
  stemr_test<lapack_complex_double, double >();
}

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
void heevd_testall_variants(){
  heevd_test<lapack_complex_float, float>();
  heevd_test<lapack_complex_double, double >();
}
void syevr_testall_variants(){
  syevr_test<float>();
  syevr_test<double>();
}
void heevr_testall_variants(){
  heevr_test<lapack_complex_float, float>();
  heevr_test<lapack_complex_double, double >();
}
void bdsqr_testall_variants(){
  bdsqr_test<float>();
  bdsqr_test<double>();
  bdsqr_test<lapack_complex_float, float>();
  bdsqr_test<lapack_complex_double, double >();
}
void bdsdc_testall_variants(){
  bdsdc_test<float>();
  bdsdc_test<double>();
}
void gesvd_testall_variants(){
  gesvd_test<float>();
  gesvd_test<float>();
  gesvd_test<lapack_complex_float, float>();
  gesvd_test<lapack_complex_double, double>();
}
void gesdd_testall_variants(){
  gesdd_test<float>();
  gesdd_test<float>();
  gesdd_test<lapack_complex_float, float>();
  gesdd_test<lapack_complex_double, double>();
}
void laswp_testall_variants(){
  laswp_test<float>();
  laswp_test<float>();
  laswp_test<lapack_complex_float>();
  laswp_test<lapack_complex_double>();
}

void laset_testall_variants(){
  laset_test<float>();
  laset_test<float>();
  laset_test<lapack_complex_float>();
  laset_test<lapack_complex_double>();
}



//#define Test(fnName)\
// test_ ## fnName ## (float)();
//
 //test_ ## fnName <double>();\
 //test_ ## fnName <lapack_complex_float>();\
 //test_ ## fnName <lapack_complex_double>();

int main(int argc, char *argv[])
{
  potrf_testall_variants();
  potf2_testall_variants();
  getrf_testall_variants(); //pivot mismatch
  getf2_testall_variants();//pivot mismatch
  geqrf_testall_variants();
  geqr2_testall_variants();
  geqpf_testall_variants();
  geqp3_testall_variants();
  gelqf_testall_variants();
  gelq2_testall_variants();
  gelsd_testall_variants();
  gelss_testall_variants();
  lauum_testall_variants();
  lauu2_testall_variants();
  potri_testall_variants();
  trtri_testall_variants();
  trti2_testall_variants();
  trsyl_testall_variants();
  gehrd_testall_variants();
  gehd2_testall_variants();
  sytrd_testall_variants();
  hetrd_testall_variants();
  sytd2_testall_variants();
  hetd2_testall_variants();
  gebrd_testall_variants();
  gebd2_testall_variants();
  sygst_testall_variants();
  hegst_testall_variants();
  sygs2_testall_variants();
  hegs2_testall_variants();
  larft_testall_variants();
  larfg_testall_variants();
  larfgp_testall_variants();
  orgqr_testall_variants();
  ungqr_testall_variants();
  ormqr_testall_variants();
  unmqr_testall_variants();
  orm2r_testall_variants();
  unm2r_testall_variants();
  orglq_testall_variants();
  unglq_testall_variants();
  ormlq_testall_variants();
  unmlq_testall_variants();
  orml2_testall_variants();
  unml2_testall_variants();
  orgtr_testall_variants();
  ungtr_testall_variants();
  ormtr_testall_variants();
  unmtr_testall_variants();
  orgbr_testall_variants();
  ungbr_testall_variants();
  ormbr_testall_variants();
  unmbr_testall_variants();
  steqr_testall_variants();
  stedc_testall_variants();
  stemr_testall_variants();
  syev_testall_variants();
  heev_testall_variants();
  syevd_testall_variants();
  heevd_testall_variants();
  syevr_testall_variants();
  heevr_testall_variants();
  bdsqr_testall_variants();
  bdsdc_testall_variants();
  gesvd_testall_variants();
  gesdd_testall_variants();
  laswp_testall_variants(); //pass //pivot mismatch
  laset_testall_variants();
}