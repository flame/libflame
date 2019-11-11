
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

/*! @file liblame_geqpf.cc
 *  libflame_geqpf.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"
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
  int n = 512;
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


void geqpf_testall_variants(){
  geqpf_test<float>();
  geqpf_test<double>();
  geqpf_test<lapack_complex_float, float>();
  geqpf_test<lapack_complex_double, double>();
}

int main(int argc, char *argv[])
{
  geqpf_testall_variants();
}