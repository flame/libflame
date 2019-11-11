
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

/*! @file liblame_unmtr.cc
 *  libflame_unmtr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename Ta, typename Tb >
void unmtr_test()
{
  int m = 8;
  int n = 8;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  char sideCPP = 'R';
  char uplo = 'U';
  char transCPP = 'N';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  tauCPPOBuff =  new Ta [m-1];
  tauCOBuff =  new Ta [m-1];
  Tb *d =  new Tb [m];
  Tb *e =  new Tb [m-1];

    //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::unmtr( LAPACK_COL_MAJOR, &sideCPP, &uplo, &transCPP, &m, &n, aCPPIOBuff, &m, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj );

  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_apply_Q_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, aCIOObj, tauCOObj, cCIOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "unmtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete cCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &cCIOObj );
  FLA_Obj_free( &tauCOObj );
}

void unmtr_testall_variants(){
  unmtr_test<lapack_complex_float, float>();
  unmtr_test<lapack_complex_double, double >();
}

int main(int argc, char *argv[])
{
  unmtr_testall_variants();
}