
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

/*! @file liblame_ormtr.cc
 *  libflame_ormtr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void ormtr_test()
{
  int m = 512;
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  char sideCPP = 'L';
  char uplo = 'L';
  char transCPP = 'T';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [m-1];
  tauCOBuff =  new T [m-1];
  T *d =  new T [m];
  T *e =  new T [m-1];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::ormtr( LAPACK_COL_MAJOR, &sideCPP, &uplo, &transCPP, &m, &n, aCPPIOBuff, &m, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_apply_Q_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormtr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormtr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] tauCPPOBuff ;
  delete[] d ;
  delete[] e ;
  delete[] aCIOBuff ;
  delete[] tauCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauCOObj );
}

void ormtr_testall_variants(){
  ormtr_test<float>();
  ormtr_test<double>();
}

int main(int argc, char *argv[])
{
  ormtr_testall_variants();
}