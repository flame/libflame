
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

/*! @file liblame_ormqr.cc
 *  libflame_ormqr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void ormqr_test()
{
  int m = 64;//512;
  int n = 128;//256;
  int k = 64;//128;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *cCPPIOBuff, *cCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();
  char side = 'R';
  char trans = 'T';

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*k);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, aCPPIOBuff, &n, tauCPPOBuff );
  libflame::ormqr( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &n, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, k, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( FLA_RIGHT, FLA_TRANSPOSE, FLA_COLUMNWISE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<T>( n, k, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormqr(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "ormqr(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] tauCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] tauCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauCOObj );

}

void ormqr_testall_variants(){
  ormqr_test<float>();
  ormqr_test<double>();
}

int main(int argc, char *argv[])
{
  ormqr_testall_variants();
}