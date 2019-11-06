
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

/*! @file liblame_gesdd.cc
 *  libflame_gesdd.cc Test application to validate CPP template interface
 *  */
#include "libflame_test.hh"

template< typename T >
void gesdd_test()
{
  int m = 64;
  int n = 64;
  int min_m_n = min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  T *aCPPIOBuff, *aCIOBuff;
  T *sCPPOBuff, *sCOBuff;
  T *uCPPOBuff, *uCOBuff;
  T *vtCPPOBuff, *vtCOBuff;

  int datatype = getDatatype<T>();
  char jobz = 'O'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new T[min_m_n];
  sCOBuff = new T[min_m_n];
  uCPPOBuff = new T[m*m];
  uCOBuff = new T[m*m];
  vtCPPOBuff = new T[m*m];
  vtCOBuff = new T[m*m];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //int gesdd( int matrix_layout, char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt )
  //Call CPP function
  libflame::gesdd( LAPACK_COL_MAJOR, &jobz, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatype, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svdd_external( FLA_SVD_VECTORS_MIN_OVERWRITE, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<T>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<T>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gessd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gessd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}

template< typename Ta, typename Tb >
void gesdd_test()
{
  int m = 64;
  int n = 64;
  int min_m_n = min(m, n);
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, sCOObj, uCOObj, vtCOObj;
  Ta *aCPPIOBuff, *aCIOBuff;
  Tb *sCPPOBuff, *sCOBuff;
  Ta *uCPPOBuff, *uCOBuff;
  Ta *vtCPPOBuff, *vtCOBuff;

  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  char jobz = 'O'; //A , S, O, N
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);

  sCPPOBuff = new Tb[min_m_n];
  sCOBuff = new Tb[min_m_n];
  uCPPOBuff = new Ta[m*m];
  uCOBuff = new Ta[m*m];
  vtCPPOBuff = new Ta[m*m];
  vtCOBuff = new Ta[m*m];
  for(int i =0; i< m*m; i++)
  {
     uCPPOBuff[i] = 0;
     uCOBuff[i] = 0;
  }

  //int gesdd( int matrix_layout, char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt )
  //Call CPP function
  libflame::gesdd( LAPACK_COL_MAJOR, &jobz, &m, &n, aCPPIOBuff, &m, sCPPOBuff, uCPPOBuff, &m, vtCPPOBuff, &n );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, min_m_n, 1, &sCOObj );
  FLA_Obj_create_without_buffer( datatype, m, m, &uCOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &vtCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( uCOBuff, 1, m, &uCOObj );
  FLA_Obj_attach_buffer( vtCOBuff, 1, n, &vtCOObj );

  //Call C function
  FLA_Svdd_external( FLA_SVD_VECTORS_MIN_OVERWRITE, aCIOObj, sCOObj, uCOObj, vtCOObj );

  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( min_m_n, 1, sCOBuff, sCPPOBuff );
  diff +=  computeError<Ta>( m, m, uCOBuff, uCPPOBuff );
  diff +=  computeError<Ta>( n, n, vtCOBuff, vtCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gessd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gessd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete sCPPOBuff ;
  delete uCPPOBuff ;
  delete vtCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &sCOObj );
  FLA_Obj_free( &uCOObj );
  FLA_Obj_free( &vtCOObj );
}

void gesdd_testall_variants(){
  gesdd_test<float>();
  gesdd_test<float>();
  gesdd_test<lapack_complex_float, float>();
  gesdd_test<lapack_complex_double, double>();
}

int main(int argc, char *argv[])
{
  gesdd_testall_variants();
}