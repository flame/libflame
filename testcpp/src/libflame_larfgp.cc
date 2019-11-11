
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

/*! @file liblame_largp.cc
 *  libflame_largp.cc Test application to validate CPP template interfaces
 *  for all libflame modules
 *  */
 
#include "libflame_test.hh"

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
    printf( "larfgp(): Failure Diff = %E\n", diff);
  }
  else{
    printf( "larfgp(): Success\n");
  }

  //Free up the buffers
  delete xCPPIOBuff ;
  FLA_Obj_free( &xCIOObj );
}

void larfgp_testall_variants(){
  larfgp_test<float>();
  larfgp_test<double>();
  larfgp_test<lapack_complex_float>();
  larfgp_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  larfgp_testall_variants();
}