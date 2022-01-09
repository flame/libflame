
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

/*! @file liblame_bdsdc.cc
 *  libflame_bdsdc.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void bdsdc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, uCOObj, vtCIObj;
  T *dCPPIOBuff, *dCIOBuff;
  T *eCPPIOBuff, *eCIOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPIBuff, *vtCIBuff;
  T *qCPPOBuff;
  int *iqCPPOBuff;
  int SMLSIZ =25;
  int LDQ = n*(11 + 2*SMLSIZ + 8*(log2(n/(SMLSIZ+1))));
  int  LDIQ = n*(3 + 3*(log2(n/(SMLSIZ+1))));
  int datatype = getDatatype<T>();
  char uplo = 'L';
  char compq = 'I'; //N, P, I
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(vtCPPIBuff, vtCIBuff, n*n);

  uCPPOBuff = new T[n*n];
  uCOBuff = new T[n*n];
  qCPPOBuff = new T[LDQ];
  iqCPPOBuff = new int[LDIQ];

  //Call CPP function
  libflame::bdsdc( LAPACK_COL_MAJOR, &uplo, &compq, &n, dCPPIOBuff, eCPPIOBuff, uCPPOBuff, &n, vtCPPIBuff, &n, qCPPOBuff, iqCPPOBuff);

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
  FLA_Bsvdd_external( FLA_LOWER_TRIANGULAR, dCIOObj, eCIOObj, uCOObj, vtCIObj );

  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCIBuff, vtCPPIBuff );

  if(diff != 0.0)
  {
    printf( "bdsdc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "bdsdc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] dCPPIOBuff ;
  delete[] eCPPIOBuff ;
  delete[] uCPPOBuff ;
  delete[] vtCPPIBuff ;
  delete[] qCPPOBuff ;
  delete[] iqCPPOBuff ;
  delete[] dCIOBuff ;
  delete[] eCIOBuff ;
  delete[] uCOBuff ;
  delete[] vtCIBuff ; 
  FLA_Obj_free_without_buffer( &dCIOObj );
  FLA_Obj_free_without_buffer( &eCIOObj );
  FLA_Obj_free_without_buffer( &uCOObj );
  FLA_Obj_free_without_buffer( &vtCIObj );
}

void bdsdc_testall_variants(){
  bdsdc_test<float>();
  bdsdc_test<double>();
}

int main(int argc, char *argv[])
{
  bdsdc_testall_variants();
}