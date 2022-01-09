
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

/*! @file liblame_getrf.cc
 *  libflame_getrd.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"

template< typename T >
void getrf_test()
{

  int m = 512;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, pivCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int *pivCPPOBuff, *pivCOBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  pivCPPOBuff =  new int [min_m_n];
  pivCOBuff =  new int [min_m_n];

 //Call CPP function
  libflame::getrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, pivCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( pivCOBuff, 1, min_m_n, &pivCOObj );

  //Call C function
  FLA_LU_piv_blk_external( aCIOObj, pivCOObj );
  FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, pivCOObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivCOBuff, pivCPPOBuff );

  if(diff != 0.0)
  {
    printf( "getrf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else if(diffInt !=0){
    printf( "getrf(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
    printf( "getrf(): Failure Pivot Buffer mismatach Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] pivCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] pivCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &pivCOObj );
}

void getrf_testall_variants(){
  getrf_test<float>();
  getrf_test<double>();
  getrf_test<lapack_complex_float>();
  getrf_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  getrf_testall_variants();
}
