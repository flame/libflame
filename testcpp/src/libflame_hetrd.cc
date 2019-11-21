
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

/*! @file liblame_hetrd.cc
 *  libflame_hetrd.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"

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

void hetrd_testall_variants(){
  hetrd_test<lapack_complex_float, float>();
  hetrd_test<lapack_complex_double, double>();
}

int main(int argc, char *argv[])
{
  hetrd_testall_variants();
}