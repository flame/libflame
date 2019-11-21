
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

/*! @file liblame_sytd2.cc
 *  libflame_sytd2.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

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

void sytd2_testall_variants(){
  sytd2_test<float>();
  sytd2_test<double>();
}

int main(int argc, char *argv[])
{
  sytd2_testall_variants();
}