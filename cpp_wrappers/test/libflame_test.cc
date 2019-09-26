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

/*! @file liblame_test.cc
 *  libflame_test.cc Test application to validate CPP template interfaces
 *  for all libflame modules
 *  */

#include "libflame_test.hh"


template< typename T >
void getrf_test()
{
  int rValCPP;
  int m = 1000;
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRrefObj, pivRefObj;
  T *aInBuff, *aRefBuff ;
  int *pivBuff, *pivRefBuff ;
  int min_m_n = min( m, n );

  int dataType = getDataType<T>();

  //Allocate and initialize the buffer with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n);
  allocate_init_buffer(pivBuff, pivRefBuff, min_m_n);

  //Allocate Object for C function and copy the buffer
  FLA_Obj_create_without_buffer( dataType, m, n, &aRrefObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, min_m_n, &aRrefObj );
  FLA_Obj_attach_buffer( pivRefBuff, 1, min_m_n, &pivRefObj );

  //Call CPP function
  libflame::getrf( &m, &n, aInBuff, &min_m_n, pivBuff, &rValCPP );

  //Call C function
  FLA_LU_piv_blk_external( aRrefObj, pivRefObj );

  //Compute Difference between C and CPP buffer
  double diff =  computeError<T>( n, m, aRefBuff, aInBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivRefBuff, pivBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getrf(): Failure Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete pivBuff ;
  FLA_Obj_free( &aRrefObj );
  FLA_Obj_free( &pivRefObj );
}

template< typename T >
void potrf_test()
{
  int rValCPP;
  int m = 1000;
  char blas_uplo;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRrefObj, pivRefObj;
  T *aInBuff, *aRefBuff ;
  int *pivBuff, *pivRefBuff ;

  int dataType = getDataType<T>();

  //Allocate and initialize the buffer with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);
  allocate_init_buffer(pivBuff, pivRefBuff, m);

  //Allocate Object for C function and copy the buffer
  FLA_Obj_create_without_buffer( dataType, m, m, &aRrefObj );
  FLA_Obj_create_without_buffer( FLA_INT, m, 1, &pivRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRrefObj );
  FLA_Obj_attach_buffer( pivRefBuff, 1, m, &pivRefObj );
  FLA_Param_map_flame_to_netlib_uplo( FLA_LOWER_TRIANGULAR, &blas_uplo );

  //Call CPP function
  libflame::potrf(&blas_uplo, &m, aInBuff, &m, &rValCPP );

  //Call C function
  FLA_Chol_blk_external( FLA_LOWER_TRIANGULAR, aRrefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );
  int diffInt =  computeError<int>( 1, m, pivRefBuff, pivBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "potrf(): Failure Diff = %d\n", diffInt);
  }else{
    printf( "potrf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  delete pivBuff;
  FLA_Obj_free( &aRrefObj );
  FLA_Obj_free( &pivRefObj );
}

void getrf_testAll(){
  getrf_test<float>();
  getrf_test<double>();
  getrf_test<scomplex>();
  getrf_test<dcomplex>();
}
void potrf_testAll(){
  potrf_test<float>();
  potrf_test<double>();
  potrf_test<scomplex>();
  potrf_test<dcomplex>();
}

int main(int argc, char *argv[])
{
  getrf_testAll();
  potrf_testAll();
}
