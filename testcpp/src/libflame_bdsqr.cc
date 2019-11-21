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

/*! @file liblame_bdsqr.cc
 *  libflame_bdsqr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void bdsqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  T *dCPPIOBuff, *dCIOBuff;
  T *eCPPIOBuff, *eCIOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPIBuff, *vtCIBuff;
  T *cCPPOBuff;
  int datatype = getDatatype<T>();
  char uplo = 'L';
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new T[n*n];
  uCOBuff = new T[n*n];
  cCPPOBuff = new T[n*n];
  int ncvt = n, nru = n, ncc = n, ldc = n;
  for(int i =0; i < n*n ; i++)
  {
    uCOBuff[i] = 0;
    uCPPOBuff[i] = 0;
  }


  //Call CPP function
  libflame::bdsqr( LAPACK_COL_MAJOR, &uplo, &n, &ncvt, &nru, &ncc, dCPPIOBuff, eCPPIOBuff, vtCPPIBuff, &n,  uCPPOBuff, &n, cCPPOBuff, &ldc);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCIObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, n, &uCOObj );
  FLA_Obj_attach_buffer( vtCIBuff, 1, n, &vtCIObj );

  //Call C function
  FLA_Bsvd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete uCPPOBuff ;
  delete vtCPPIBuff ;
  delete cCPPOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCIObj );
}

template< typename Ta, typename Tb >
void bdsqr_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  Tb *dCPPIOBuff, *dCIOBuff;
  Tb *eCPPIOBuff, *eCIOBuff;
  Ta *uCPPOBuff, *uCOBuff;
  Ta *vtCPPIBuff, *vtCIBuff;
  Ta *cCPPOBuff;
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  char uplo = 'L';
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new Ta[n*n];
  uCOBuff = new Ta[n*n];
  cCPPOBuff = new Ta[n*n];
  int ncvt = n, nru = n, ncc = n, ldc = n;
  //Call CPP function
  libflame::bdsqr( LAPACK_COL_MAJOR, &uplo, &n, &ncvt, &nru, &ncc, dCPPIOBuff, eCPPIOBuff, vtCPPIBuff, &n,  uCPPOBuff, &n, cCPPOBuff, &ldc);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &uCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, n, &vtCIObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, n, &uCOObj );
  FLA_Obj_attach_buffer( vtCIBuff, 1, n, &vtCIObj );

  //Call C function
  FLA_Bsvd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<Tb>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<Tb>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<Ta>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<Ta>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete uCPPOBuff ;
  delete vtCPPIBuff ;
  delete cCPPOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCIObj );
}

void bdsqr_testall_variants(){
  bdsqr_test<float>();
  bdsqr_test<double>();
  bdsqr_test<lapack_complex_float, float>();
  bdsqr_test<lapack_complex_double, double >();
}

int main(int argc, char *argv[])
{
  bdsqr_testall_variants();
}