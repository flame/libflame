
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

/*! @file liblame_hegs2.cc
 *  libflame_hegs2.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void hegs2_test()
{
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 2 ; //1 or 2
  int itype_c = FLA_NO_INVERSE;

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);


  //Call CPP function
  libflame::hegs2( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bCIOObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bCIOObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegs2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "hegs2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
}

void hegs2_testall_variants(){
  hegs2_test<lapack_complex_float>();
  hegs2_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  hegs2_testall_variants();
}