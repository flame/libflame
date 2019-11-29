
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

/*! @file liblame_heevr.cc
 *  libflame_heevr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename Ta, typename Tb >
void heevr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj;
  Ta *aCPPIOBuff, *zCPPOBuff;
  Tb  *wCPPOBuff;
  int *isuppzCPPOBuff;
  Ta *aCIOBuff, *zCOBuff;
  Tb *wCOBuff;
  char jobz = 'N'; //N, V
  char uplo = 'l';
  char range = 'A'; //A, V, I
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  Tb v1 = 0, vu = 0, abstol = FLA_MACH_SFMIN;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];
  zCOBuff  = new Ta[n*n];
  zCPPOBuff = new Ta[n*n];
  isuppzCPPOBuff = new int [isuppzDim];

  for(int i =0; i < n*n ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::heevr( LAPACK_COL_MAJOR, &jobz, &range, &uplo, &n, aCPPIOBuff, &n, &v1, &vu, &il, &iu,
                   &abstol, &m, wCPPOBuff, zCPPOBuff, &ldz, isuppzCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &zCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );

  //Call C function
  FLA_Hevdr_external( FLA_EVD_WITHOUT_VECTORS, FLA_LOWER_TRIANGULAR, aCIOObj, wCOObj, zCOObj );

  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<Ta>( n, n, zCOBuff, zCPPOBuff );

  if(diff != 0.0)
  {
    printf( "heevr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "heevr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] wCPPOBuff ;
  delete[] zCPPOBuff ;
  delete[] isuppzCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] wCOBuff ;
  delete[] zCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &wCOObj );
  FLA_Obj_free_without_buffer( &zCOObj );
}

void heevr_testall_variants(){
  heevr_test<lapack_complex_float, float>();
  heevr_test<lapack_complex_double, double >();
}

int main(int argc, char *argv[])
{
  heevr_testall_variants();
}