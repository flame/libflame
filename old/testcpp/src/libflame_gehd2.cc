
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

/*! @file liblame_gehd2.cc
 *  libflame_gehd2.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

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
  delete[] aCPPIOBuff ;
  delete[] tauCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] tauCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauCOObj );
}

void gehd2_testall_variants(){
  gehd2_test<float>();
  gehd2_test<double>();
  gehd2_test<lapack_complex_float>();
  gehd2_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  gehd2_testall_variants();
}