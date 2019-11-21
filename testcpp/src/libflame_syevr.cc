
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

/*! @file liblame_syevr.cc
 *  libflame_syevr.cc Test application to validate CPP template interface
 *  */

#include "libflame_test.hh"

template< typename T >
void syevr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj;
  T *aCPPIOBuff, *wCPPOBuff, *zCPPOBuff;
  int *isuppzCPPOBuff;
  T *aCIOBuff, *wCOBuff, *zCOBuff;
  char jobz = 'N'; //N, V
  char uplo = 'l';
  char range = 'A'; //A, V, I
  int datatype = getDatatype<T>();
  T v1 = 0, vu = 0, abstol = FLA_MACH_SFMIN;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  zCOBuff  = new T[n*n];
  zCPPOBuff = new T[n*n];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*n ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::syevr( LAPACK_COL_MAJOR, &jobz, &range, &uplo, &n, aCPPIOBuff, &n, &v1, &vu, &il, &iu,
                   &abstol, &m, wCPPOBuff, zCPPOBuff, &ldz, isuppzCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );

  //Call C function
  FLA_Hevdr_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, aCIOObj, wCOObj, zCOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<T>( n, n, zCOBuff, zCPPOBuff );

  if(diff != 0.0)
  {
    printf( "syevr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "syevr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
}

void syevr_testall_variants(){
  syevr_test<float>();
  syevr_test<double>();
}

int main(int argc, char *argv[])
{
  syevr_testall_variants();
}
