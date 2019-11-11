
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

/*! @file liblame_stemr.cc
 *  libflame_stemr.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void stemr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj, eCIOObj;
  T *aCPPIOBuff, *wCPPOBuff, *zCPPOBuff, *eCPPIOBuff;
  int *isuppzCPPOBuff;
  T *aCIOBuff, *wCOBuff, *zCOBuff, *eCIOBuff;
  char jobz = 'N'; //N, V
  char range = 'A'; //A, V, I
  int datatype = getDatatype<T>();
  int v1 = 0, vu = 0;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  int tryrac = 1;
  int nzc = n;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, n);
  wCPPOBuff  = new T[n];
  wCOBuff  = new T[n];
  zCOBuff  = new T[n*m];
  zCPPOBuff = new T[n*m];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*m ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::stemr( LAPACK_COL_MAJOR, &jobz, &range, &n, aCPPIOBuff, eCPPIOBuff, &v1, &vu, &il, &iu,
                   &m, wCPPOBuff, zCPPOBuff, &ldz, &nzc, isuppzCPPOBuff, &tryrac);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatype, n, m, &zCOObj );
  FLA_Obj_create_without_buffer( datatype, n, 1, &eCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, n, &eCIOObj );

  //Call C function
  FLA_Tevdr_external( FLA_EVD_WITHOUT_VECTORS, aCIOObj, eCIOObj, wCOObj, zCOObj );

  double diff =  computeError<T>( n, 1, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<T>( n, m, zCOBuff, zCPPOBuff );
  diff +=  computeError<T>( n, 1, eCIOBuff, eCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stemr(): Failure Diff = %E\n", diff);
  }else{
    printf( "stemr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete eCPPIOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
  FLA_Obj_free( &eCIOObj );
}

template< typename Ta, typename Tb >
void stemr_test()
{
  int n = 512;
  int m = n;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, wCOObj, zCOObj, eCIOObj;
   Ta *zCPPOBuff ;
  Tb *wCPPOBuff, *aCPPIOBuff, *eCPPIOBuff;
  int *isuppzCPPOBuff;
  Ta *zCOBuff;
  Tb *wCOBuff, *aCIOBuff, *eCIOBuff;
  char jobz = 'N'; //N, V
  char range = 'A'; //A, V, I
  int datatypeTa = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  int v1 = 0, vu = 0;
  int il = 0 , iu = 0 ;
  int ldz = 1; //1, n
  int isuppzDim = 2*max(1,m) ;
  int tryrac = 1;
  int nzc = n;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, n);
  wCPPOBuff  = new Tb[n];
  wCOBuff  = new Tb[n];
  zCOBuff  = new Ta[n*m];
  zCPPOBuff = new Ta[n*m];
  isuppzCPPOBuff = new int [isuppzDim];

   for(int i =0; i < n*m ; i++)
  {
    zCPPOBuff[i] = 0;
    zCOBuff[i] = 0;
  }

  //Call CPP function
  libflame::stemr( LAPACK_COL_MAJOR, &jobz, &range, &n, aCPPIOBuff, eCPPIOBuff, &v1, &vu, &il, &iu,
                   &m, wCPPOBuff, zCPPOBuff, &ldz, &nzc, isuppzCPPOBuff, &tryrac);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &aCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &wCOObj );
  FLA_Obj_create_without_buffer( datatypeTa, n, m, &zCOObj );
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &eCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( wCOBuff, 1, n, &wCOObj );
  FLA_Obj_attach_buffer( zCOBuff, 1, n, &zCOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, n, &eCIOObj );

  //Call C function
  FLA_Tevdr_external( FLA_EVD_WITHOUT_VECTORS, aCIOObj, eCIOObj, wCOObj, zCOObj );

  double diff =  computeError<Tb>( n, 1, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Tb>( n, 1, wCOBuff, wCPPOBuff );
  diff +=  computeError<Ta>( n, m, zCOBuff, zCPPOBuff );
  diff +=  computeError<Tb>( n, 1, eCIOBuff, eCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stemr(): Failure Diff = %E\n", diff);
  }else{
    printf( "stemr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete wCPPOBuff ;
  delete zCPPOBuff ;
  delete eCPPIOBuff ;
  delete isuppzCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &wCOObj );
  FLA_Obj_free( &zCOObj );
  FLA_Obj_free( &eCIOObj );
}

void stemr_testall_variants(){
  stemr_test<float>();
  stemr_test<double>();
  stemr_test<lapack_complex_float, float>();
  stemr_test<lapack_complex_double, double >();
}

int main(int argc, char *argv[])
{
  stemr_testall_variants();
}