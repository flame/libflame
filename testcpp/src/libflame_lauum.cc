
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

/*! @file liblame_lauum.cc
 *  libflame_lauum.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

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

void lauum_testall_variants(){
  lauum_test<float>();
  lauum_test<double>();
  lauum_test<lapack_complex_float>();
  lauum_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  lauum_testall_variants();
}