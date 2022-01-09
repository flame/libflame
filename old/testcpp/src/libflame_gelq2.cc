
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

/*! @file liblame_gelq2.cc
 *  libflame_gelq2.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"

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
  int min_m_n = min( m, n ) ;
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
  delete[] aCPPIOBuff ;
  delete[] tauCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] tauCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauCOObj );
}

void gelq2_testall_variants(){
  gelq2_test<float>();
  gelq2_test<double>();
  gelq2_test<lapack_complex_float>();
  gelq2_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  gelq2_testall_variants();
}