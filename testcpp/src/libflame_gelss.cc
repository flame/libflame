
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

/*! @file liblame_gelss.cc
 *  libflame_gelss.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"

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


void gelss_testall_variants(){
  gelss_test<float>();
  gelss_test<double>();
  gelss_test<lapack_complex_float, float>();
  gelss_test<lapack_complex_double, double>();
}

int main(int argc, char *argv[])
{
  gelss_testall_variants();
}