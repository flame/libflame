
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

/*! @file liblame_unglq.cc
 *  libflame_unglq.cc Test application to validate CPP template interface
 *  */
#include "libflame_test.hh"

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
  int min_m_n = min(m,n);

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
    printf( "unglq(): Failure Diff = %E\n", diff);
  }else{
    printf( "unglq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

void unglq_testall_variants(){
  unglq_test<lapack_complex_float>();
  unglq_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  unglq_testall_variants();
}