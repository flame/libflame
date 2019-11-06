
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

/*! @file liblame_ungqr.cc
 *  libflame_ungqr.cc Test application to validate CPP template interface
 *  */
#include "libflame_test.hh"

template< typename T >
void ungqr_test()
{
  int m = 512;
  int n = 256;
  int k = 256;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
  libflame::ungqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj );

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_QR_form_Q_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "ungqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

void ungqr_testall_variants(){
  ungqr_test<lapack_complex_float>();
  ungqr_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  ungqr_testall_variants();
}