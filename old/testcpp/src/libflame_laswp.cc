
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

/*! @file liblame_laswp.cc
 *  libflame_laswp.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void laswp_test()
{
  int n = 1024;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, ipivIObj;
  T *aCPPIOBuff, *aCIOBuff;
  int *ipivCPPIBuff, *ipivCIBuff;
  int k1 = 1;
  int k2 =  128;
  int incx = 1;
  int pivDim = (k1+(k2-k1)*abs(incx));
  int datatype = getDatatype<T>();
  int datatypeTb = getDatatype<int>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(ipivCPPIBuff, ipivCIBuff, pivDim);

  for (int i = 0; i < pivDim; i++ )
  {
    ipivCPPIBuff[ i ] = ipivCIBuff[ i ] + i + 1;
  }

  //Call CPP function
  libflame::laswp( LAPACK_COL_MAJOR, &n, aCPPIOBuff, &n, &k1, &pivDim, ipivCPPIBuff, &incx );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, pivDim, 1, &ipivIObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( ipivCIBuff, 1, pivDim, &ipivIObj );

  //Call C function
  FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, ipivIObj, aCIOObj  );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "laswp(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "laswp(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] ipivCPPIBuff ;
  delete[] aCIOBuff ;
  delete[] ipivCIBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &ipivIObj );
}

void laswp_testall_variants(){
  laswp_test<float>();
  laswp_test<float>();
  laswp_test<lapack_complex_float>();
  laswp_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  laswp_testall_variants();
}