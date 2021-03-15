
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

/*! @file liblame_laset.cc
 *  libflame_laset.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

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
  delete[] aCPPIOBuff ;
  delete[] aCIOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
}

void laset_testall_variants(){
  laset_test<float>();
  laset_test<float>();
  laset_test<lapack_complex_float>();
  laset_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  laset_testall_variants();
}