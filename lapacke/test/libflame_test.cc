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
typedef int FLA_Error;

template< typename T >
void potrf_test()
{

  int m = 32 ;
  char blas_uplo = 'U';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::potrf(LAPACK_COL_MAJOR, &blas_uplo, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );  

  //Call C function
  FLA_Chol_blk_external( FLA_UPPER_TRIANGULAR, aRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "potrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "potrf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}

template< typename T >
void potf2_test()
{

  int m = 16; //2048
  char blas_uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::potf2( LAPACK_COL_MAJOR, &blas_uplo, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  
  //Call C function
  FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, aRefObj );
 
  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "potf2():Failure Diff = %E\n", diff);
  }else{
    printf( "potf2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}

template< typename T >
void getrf_test()
{

  int m = 512;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, pivRefObj;
  T *aInBuff, *aRefBuff ;
  int *pivBuff, *pivRefBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n);
  pivBuff =  new int [min_m_n];
  pivRefBuff =  new int [min_m_n];


 //Call CPP function
  libflame::getrf( LAPACK_COL_MAJOR, &m, &n, aInBuff, &m, pivBuff );
 
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( pivRefBuff, 1, min_m_n, &pivRefObj );
  //FLA_Set( FLA_ZERO, pivRefObj );  
  
  //Call C function
  FLA_LU_piv_blk_external( aRefObj, pivRefObj );


  double diff =  computeError<T>( n, m, aRefBuff, aInBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivRefBuff, pivBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getrf(): Failure Pivot BUffer mismatach Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete pivBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &pivRefObj );
}


template< typename T >
void getf2_test()
{

  int m = 1024;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, pivRefObj;
  T *aInBuff, *aRefBuff ;
  int *pivBuff, *pivRefBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n);
  pivBuff =  new int [min_m_n];
  pivRefBuff =  new int [min_m_n];

 //Call CPP function
  libflame::getf2( LAPACK_COL_MAJOR, &m, &n, aInBuff, &m, pivBuff );
 
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( pivRefBuff, 1, min_m_n, &pivRefObj );
  //FLA_Set( FLA_ZERO, pivRefObj );  
  
  //Call C function
  FLA_LU_piv_unb_external( aRefObj, pivRefObj );

  double diff =  computeError<T>( n, m, aRefBuff, aInBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivRefBuff, pivBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getf2(): Failure:Pivot Buffer Mismatch Diff = %d\n", diffInt);
  }else{
    printf( "getf2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete pivBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &pivRefObj );
}
#if 0
#define PREFIX2FLAME_INVERT_TAU(datatype, tauRefObj)\
{\
  if (datatype == FLA_FLOAT)\
    FLAME_invert_stau(tauRefObj);\
  if (datatype == FLA_DOUBLE)\
    FLAME_invert_dtau(tauRefObj);\
  if (datatype == FLA_COMPLEX)\
    FLAME_invert_ctau(tauRefObj);\
  if (datatype == FLA_DOUBLE_COMPLEX)\
    FLAME_invert_ztau(tauRefObj);\
}
#endif

template< typename T >
void geqrf_test()
{

  int rValC = 0;
  int m = 512;
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, tObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int min_m_n = min( m, n );
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n);
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

 //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ); 

  //Call C function
  rValC = FLA_QR_blk_external( aRefObj, tauRefObj );


   fp = fopen("test/out","a+");
   print(fp, m*n,aInBuff, aRefBuff);
   fclose(fp);
  
  double diff = computeError<T>( n, m, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff );

  if(diff != 0.0)
  {
    printf( "geqrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "geqrf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void geqr2_test()
{

  int rValC;
  int m = 128;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [m];

  //Call CPP function
  libflame::geqr2( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, tauBuff );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;
 
  //Call C function
  rValC = FLA_QR_unb_external( aRefObj, tauRefObj );



  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqr2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqr2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}
#if 0
template< typename T >
void geqpf_test()
{

  int rValC;
  int m = 256;
  int n = 4096;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int *jpvtBuff, *jpvtRefBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [3*n];
  jpvtBuff = new int[n];
  jpvtRefBuff = new int[n];
  
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;

  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, jpvtBuff, tauBuff );

  //Call C function
  //rValC = geqpf_C( aRefObj, jpvtRefBuff, tauRefObj, rworkObj );

  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename Ta, typename Tb >
void geqpf_test()
{

  int rValC;
  int m = 512;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauBuff, *tauRefBuff ;
  Ta *workBuff;
  Tb *rworkBuff, *rworkRefBuff;
  int *jpvtBuff, *jpvtRefBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  tauBuff =  new Ta [min_m_n];
  tauRefBuff =  new Ta [min_m_n];
  workBuff =  new Ta [3 * n];
  rworkBuff =  new Tb [2 * n];
  rworkRefBuff = new Tb [2 * n];
  jpvtBuff = new int[n];
  jpvtRefBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  allocate_init_buffer(rworkBuff, rworkRefBuff, 2*n,0) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;
  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, jpvtBuff, tauBuff ) ;

  //Call C function
  //rValC = geqpf_C( aRefObj, jpvtRefBuff, tauRefObj, rworkRefObj );

  
  double diff =  computeError<Ta>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<Ta>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  diff +=  computeError<Tb>( 1, 2*n, rworkRefBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "geqpf(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete rworkBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
  FLA_Obj_free( &rworkRefObj );
}



template< typename T >
void geqp3_test()
{

  int rValC;
  int m = 256;
  int n = 4096;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int *jpvtBuff, *jpvtRefBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [3*n];
  jpvtBuff = new int[n];
  jpvtRefBuff = new int[n];
  
  int datatype = getDatatype<T>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, jpvtBuff, tauBuff );

  //Call C function
  //rValC = geqp3_C( aRefObj, jpvtRefBuff, tauRefObj, rworkObj );


  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename Ta, typename Tb >
void geqp3_test()
{

  int rValC;
  int m = 2048;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauBuff, *tauRefBuff ;
  Ta *workBuff;
  Tb *rworkBuff, *rworkRefBuff;
  int *jpvtBuff, *jpvtRefBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  tauBuff =  new Ta [min_m_n];
  tauRefBuff =  new Ta [min_m_n];
  workBuff =  new Ta [3 * n];
  rworkBuff =  new Tb [2 * n];
  rworkRefBuff = new Tb [2 * n];
  jpvtBuff = new int[n];
  jpvtRefBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  allocate_init_buffer(rworkBuff, rworkRefBuff, 2*n,0) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, jpvtBuff, tauBuff ) ;

  //Call C function
  rValC = geqp3_C( aRefObj, jpvtRefBuff, tauRefObj, rworkRefObj );

  //Check for errors of C and CPP function call
  if( rValCPP !=0 || rValC != 0){
    printf("geqp3(): Status of C algo: %d CPP algo: %d\n", rValC, rValCPP);
    return;
  }
  
  double diff =  computeError<Ta>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<Ta>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  diff +=  computeError<Tb>( 1, 2*n, rworkRefBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "geqp3(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete rworkBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
  FLA_Obj_free( &rworkRefObj );
}

#endif
template< typename T >
void gelqf_test()
{

  int rValC;
  int m = 8192;
  int n = 1024;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int min_m_n = min( m, n ) ;

  int datatype = getDatatype<T>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &m, &n, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;

  //Call C function
  rValC = FLA_LQ_blk_external( aRefObj, tauRefObj );



  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelqf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelqf(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void gelq2_test()
{

  int rValC;
  int m = 1000;
  int n = 1000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [m];

  //Call CPP function
  libflame::gelq2( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;

  //Call C function
  rValC = FLA_LQ_unb_external( aRefObj, tauRefObj );
 
  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelq2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelq2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

#if 0
template< typename T >
void gelsd_test()
{

  int rValC;
  int m = 16384;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff, s, rcond ;
  int *b, *bRef,  rank;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));
  
  tauBuff =  new T [min_m_n];
  tauRefBuff =  new T [min_m_n];
  workBuff =  new T [3*n];
  b = new int[n];
  bRef = new int[n];
  
  int datatype = getDatatype<T>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  int iwork ;
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;

  //Call CPP function
  //( int* m, int* n, int* nrhs, T* a, int* lda, T* b, int* ldb, T* s, T* rcond, int* rank, T* work, int* lwork, int* iwork, int* info )
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, &n, aInBuff, &lda, b, &ldb, s, rcond, rank );

  //Call C function
  rValC = gelsd_C( aRefObj, jpvtRefBuff, tauRefObj, rworkObj );

  //Check for errors of C and CPP function call
  if( rValCPP !=0 || rValC != 0){
    printf("gelsd(): Status of C algo: %d CPP algo: %d\n", rValC, rValCPP);
    return;
  }
  double diff =  computeError<T>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelsd(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelsd(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename Ta, typename Tb >
void gelsd_test()
{

  int rValC;
  int m = 2048;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj, rworkRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauBuff, *tauRefBuff ;
  Ta *workBuff;
  Tb *rworkBuff, *rworkRefBuff;
  int *jpvtBuff, *jpvtRefBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  tauBuff =  new Ta [min_m_n];
  tauRefBuff =  new Ta [min_m_n];
  workBuff =  new Ta [3 * n];
  rworkBuff =  new Tb [2 * n];
  rworkRefBuff = new Tb [2 * n];
  jpvtBuff = new int[n];
  jpvtRefBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*n) ;
  allocate_init_buffer(rworkBuff, rworkRefBuff, 2*n,0) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauRefObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, min_m_n, &tauRefObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, aInBuff, &lda, jpvtBuff, tauBuff ) ;

  //Call C function
  rValC = gelsd_C( aRefObj, jpvtRefBuff, tauRefObj, rworkRefObj );

  //Check for errors of C and CPP function call
  if( rValCPP !=0 || rValC != 0){
    printf("geqp3(): Status of C algo: %d CPP algo: %d\n", rValC, rValCPP);
    return;
  }
  
  double diff =  computeError<Ta>( n, m, aRefBuff, aInBuff ) ;
  int diffInt =  computeError<Ta>( 1, min_m_n, tauRefBuff, tauBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtRefBuff, jpvtBuff ) ;
  diff +=  computeError<Tb>( 1, 2*n, rworkRefBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "gelsd(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelsd(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelsd(): Success\n");
  }
 
  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete rworkBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
  FLA_Obj_free( &rworkRefObj );
}
#endif


template< typename T >
void lauum_test()
{

  int m = 64;
  char blas_uplo = 'u';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::lauum( LAPACK_COL_MAJOR, &blas_uplo, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );

  //Call C function
  FLA_Ttmm_blk_external( FLA_UPPER_TRIANGULAR, aRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "lauum(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauum(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}

template< typename T >
void lauu2_test()
{

  int m = 128;
  char blas_uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::lauu2( LAPACK_COL_MAJOR, &blas_uplo, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Param_map_flame_to_netlib_uplo( FLA_LOWER_TRIANGULAR, &blas_uplo );

  //Call C function
  FLA_Ttmm_unb_external( FLA_LOWER_TRIANGULAR, aRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "lauu2(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauu2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}


template< typename T >
void potri_test()
{

  int rValC;
  int m = 64; //64 works
  char blas_uplo = 'L';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::potri( LAPACK_COL_MAJOR, &blas_uplo, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );

  //Call C function
//  rValC = potri_C( FLA_LOWER_TRIANGULAR, aRefObj );



  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );
  
  if(diff != 0.0)
  {
	fp = fopen("test/in","a+");
   print(fp,m*m,aInBuff, aRefBuff);
   fclose(fp);
  }
  
  if(diff != 0.0)
  {
    printf( "potri(): Failure Diff = %E\n", diff);
  }else{
    printf( "potri(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}



template< typename T >
void trtri_test()
{

  int m = 128;
  char blas_uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::trtri( LAPACK_COL_MAJOR, &blas_uplo, &blas_diag, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  
  //Call C function
  FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "trtri(): Failure Diff = %E\n", diff);
  }else{
    printf( "trtri(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}

template< typename T >
void trti2_test()
{

  int m = 128;
  char blas_uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj;
  T *aInBuff, *aRefBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, m*m);

  //Call CPP function
  libflame::trti2( LAPACK_COL_MAJOR, &blas_uplo, &blas_diag, &m, aInBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  
  //Call C function
  FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aRefBuff, aInBuff );

  if(diff != 0.0)
  {
    printf( "trti2(): Failure Diff = %E\n", diff);
  }else{
    printf( "trti2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff;
  FLA_Obj_free( &aRefObj );
}

template< typename T >
void gehrd_test()
{

  int rValC;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  libflame::gehrd( LAPACK_COL_MAJOR, &n, &ilo, &iho, aInBuff, &n, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  rValC = FLA_Hess_blk_external( aRefObj, tauRefObj, ilo, iho );
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );

  if(diff != 0.0)
  {
    printf( "gehrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehrd(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void gehd2_test()
{

  int rValC;
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  rValC = libflame::gehd2( LAPACK_COL_MAJOR, &n, &ilo, &iho, aInBuff, &n, tauBuff );
printf("%d", rValC);
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  FLA_Hess_unb_external( aRefObj, tauRefObj, ilo, iho );

  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );
     fp = fopen("test/out1","a+");
   print(fp,n-1,tauRefBuff, tauRefBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,n-1,tauBuff, tauBuff);
   fclose(fp);
  if(diff != 0.0)
  {
    printf( "gehd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehd2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void sytrd_test()
{

  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *d, *e ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauRefBuff =  new T [n-1];
  work = new T [lwork];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &blas_uplo, &n, aInBuff, &n, d, e, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aRefObj, tauRefObj );
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );

  if(diff != 0.0)
  {
    printf( "sytrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytrd(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename Ta, typename Tb >
void hetrd_test()
{

  int rValC;
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauBuff, *tauRefBuff ;
  Tb  *d, *e ;
  Ta *work;
  char blas_uplo = 'u';
  int datatype = getDatatype<Ta>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauRefBuff =  new Ta [n-1];
  work = new Ta [lwork];
  
  //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &blas_uplo, &n, aInBuff, &n, d, e, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  rValC = FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aRefObj, tauRefObj );
  double diff =  computeError<Ta>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<Ta>( 1, n-1, tauRefBuff, tauBuff );
  
  if(diff != 0.0)
  {
    printf( "hetrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetrd(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void sytd2_test()
{

  int rValC;
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *workBuff, *d, *e ;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauRefBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::sytd2( LAPACK_COL_MAJOR, &blas_uplo, &n, aInBuff, &n, d, e, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  rValC = FLA_Tridiag_unb_external( FLA_LOWER_TRIANGULAR, aRefObj, tauRefObj );
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );
  
  if(diff != 0.0)
  {
    printf( "sytd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytd2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename Ta,  typename Tb >
void hetd2_test()
{

  int rValC;
  int n = 128;//2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauBuff, *tauRefBuff ;
  Ta *workBuff;
  Tb *d, *e ;
  char blas_uplo = 'u';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauRefBuff =  new Ta [n-1];
  workBuff =  new Ta [n];
   fp = fopen("test/in","a+");
   print(fp,n*n,aInBuff, aRefBuff);
   fclose(fp);
  //Call CPP function
  libflame::hetd2( LAPACK_COL_MAJOR, &blas_uplo, &n, aInBuff, &n, d, e, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, n-1, &tauRefObj ); 

  //Call C function
  rValC = FLA_Tridiag_unb_external( FLA_UPPER_TRIANGULAR, aRefObj, tauRefObj );
  double diff =  computeError<Ta>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<Ta>( 1, n-1, tauRefBuff, tauBuff );

     fp = fopen("test/out","a+");
   print(fp,n*n,aInBuff, aRefBuff);
   fclose(fp);
     fp = fopen("test/out1","a+");
   print(fp,n*n,aInBuff, aInBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,n*n,aRefBuff, aRefBuff);
   fclose(fp);

  if(diff != 0.0)
  {
    printf( "hetd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetd2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  delete workBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void gebrd_test()
{

  int rValC;
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauqRefObj, taupRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauqBuff, *tauqRefBuff ;
  T *taupBuff, *taupRefBuff ;
  T *d, *e ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqBuff =  new T [min_m_n];
  taupBuff =  new T [min_m_n];
  tauqRefBuff =  new T [min_m_n];
  taupRefBuff =  new T [min_m_n];
  work = new T [lwork];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aInBuff, &n, d, e, tauqBuff, taupBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauqRefBuff, 1, min_m_n, &tauqRefObj ); 
  FLA_Obj_attach_buffer( taupRefBuff, 1, min_m_n, &taupRefObj ); 

  //Call C function
  rValC = FLA_Bidiag_blk_external( aRefObj, tauqRefObj, taupRefObj );
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqRefBuff, tauqBuff );
  diff +=  computeError<T>( 1, min_m_n, taupRefBuff, taupBuff );

  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauqBuff ;
  delete taupBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauqRefObj );
  FLA_Obj_free( &taupRefObj );
}

template< typename Ta, typename Tb >
void gebrd_test()
{

  int rValC;
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauqRefObj, taupRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauqBuff, *tauqRefBuff ;
  Ta *taupBuff, *taupRefBuff ;
  Tb *d, *e ;
  Ta *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<Ta>();
  int lwork    = (m + n) * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqBuff =  new Ta [min_m_n];
  taupBuff =  new Ta [min_m_n];
  tauqRefBuff =  new Ta [min_m_n];
  taupRefBuff =  new Ta [min_m_n];
  work = new Ta [lwork];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aInBuff, &n, d, e, tauqBuff, taupBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauqRefBuff, 1, min_m_n, &tauqRefObj ); 
  FLA_Obj_attach_buffer( taupRefBuff, 1, min_m_n, &taupRefObj ); 

  //Call C function
  rValC = FLA_Bidiag_blk_external( aRefObj, tauqRefObj, taupRefObj );
  double diff =  computeError<Ta>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqRefBuff, tauqBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupRefBuff, taupBuff );
  
  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauqBuff ;
  delete taupBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauqRefObj );
  FLA_Obj_free( &taupRefObj );
}

template< typename T >
void gebd2_test()
{

  int rValC;
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauqRefObj, taupRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauqBuff, *tauqRefBuff ;
  T *taupBuff, *taupRefBuff ;
  T *d, *e ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqBuff =  new T [min_m_n];
  taupBuff =  new T [min_m_n];
  tauqRefBuff =  new T [min_m_n];
  taupRefBuff =  new T [min_m_n];
  work = new T [lwork];

  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aInBuff, &n, d, e, tauqBuff, taupBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauqRefBuff, 1, min_m_n, &tauqRefObj ); 
  FLA_Obj_attach_buffer( taupRefBuff, 1, min_m_n, &taupRefObj ); 

  //Call C function
  rValC = FLA_Bidiag_unb_external( aRefObj, tauqRefObj, taupRefObj );
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqRefBuff, tauqBuff );
  diff +=  computeError<T>( 1, min_m_n, taupRefBuff, taupBuff );

  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauqBuff ;
  delete taupBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauqRefObj );
  FLA_Obj_free( &taupRefObj );
}


template< typename Ta, typename Tb >
void gebd2_test()
{

  int rValC;
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauqRefObj, taupRefObj;
  Ta *aInBuff, *aRefBuff ;
  Ta *tauqBuff, *tauqRefBuff ;
  Ta *taupBuff, *taupRefBuff ;
  Tb *d, *e ;
  Ta *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<Ta>();
  int lwork    = (m + n) * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqBuff =  new Ta [min_m_n];
  taupBuff =  new Ta [min_m_n];
  tauqRefBuff =  new Ta [min_m_n];
  taupRefBuff =  new Ta [min_m_n];
  work = new Ta [lwork];
  fp = fopen("test/in","a+");
   print(fp,m*n,aInBuff, aRefBuff);
   fclose(fp);
  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aInBuff, &n, d, e, tauqBuff, taupBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqRefObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauqRefBuff, 1, min_m_n, &tauqRefObj ); 
  FLA_Obj_attach_buffer( taupRefBuff, 1, min_m_n, &taupRefObj ); 

  //Call C function
  rValC = FLA_Bidiag_unb_external( aRefObj, tauqRefObj, taupRefObj );
  double diff =  computeError<Ta>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqRefBuff, tauqBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupRefBuff, taupBuff );
   fp = fopen("test/tau1","a+");
   print(fp,min_m_n,tauqRefBuff, taupRefBuff);
   fclose(fp);
   fp = fopen("test/tau2","a+");
   print(fp,min_m_n,tauqBuff, taupBuff);
   fclose(fp);
     fp = fopen("test/out","a+");
   print(fp,m*n,aInBuff, aRefBuff);
   fclose(fp);
     fp = fopen("test/out1","a+");
   print(fp,m*n,aInBuff, aInBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,m*n,aRefBuff, aRefBuff);
   fclose(fp);
   
  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauqBuff ;
  delete taupBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauqRefObj );
  FLA_Obj_free( &taupRefObj );
}

template< typename T >
void sygst_test()
{

  int rValC;
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, bRefObj;
  T *aInBuff, *aRefBuff, *bRefBuff, *b;
  char blas_uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::sygst( LAPACK_COL_MAJOR, &itype, &blas_uplo, &n, aInBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  rValC = FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aRefObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygst(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete b ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void hegst_test()
{

  int rValC;
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, bRefObj;
  T *aInBuff, *aRefBuff, *bRefBuff, *b;
  char blas_uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::hegst( LAPACK_COL_MAJOR, &itype, &blas_uplo, &n, aInBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  rValC = FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aRefObj, bRefObj );

  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegst(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegst(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete b ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void sygs2_test()
{

  int n = 64;//2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, bRefObj;
  T *aInBuff, *aRefBuff, *bRefBuff, *b;
  char blas_uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::sygs2( LAPACK_COL_MAJOR, &itype, &blas_uplo, &n, aInBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aRefObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete b ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void hegs2_test()
{

  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, bRefObj;
  T *aInBuff, *aRefBuff, *bRefBuff, *b;
  char blas_uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::hegs2( LAPACK_COL_MAJOR, &itype, &blas_uplo, &n, aInBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, n, &aRefObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aRefObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegs2(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete b ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void orgqr_test()
{

  int rValC;
  int m = 128;
  int n = 256;
  int k = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  work = new T [lwork];

  //Call CPP function
  libflame::orgqr( LAPACK_COL_MAJOR, &m, &n, &k, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, k-1, &tauRefObj ); 

  //Call C function
  rValC = FLA_QR_form_Q_external( aRefObj, tauRefObj );
  double diff =  computeError<T>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );
  
  if(diff != 0.0)
  {
    printf( "orgqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgqr(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

template< typename T >
void ungqr_test()
{

  int rValC;
  int m = 512;
  int n = 256;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  work = new T [lwork];

  //Call CPP function
  libflame::ungqr( LAPACK_COL_MAJOR, &m, &n, &k, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, k-1, &tauRefObj ); 

  //Call C function
//  rValC =  ungqr_C( aRefObj, tauRefObj );
  double diff =  computeError<T>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );
  
  if(diff != 0.0)
  {
    printf( "ungqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungqr(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}

#if 0
template< typename T >
void orglq_test()
{

  int rValC;
  int m = 128;
  int n = 256;
  int k = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  work = new T [lwork];

  //Call CPP function
  libflame::orglq( &m, &n, &k, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, k-1, &tauRefObj ); 

  //Call C function
 // rValC = orglq_C( aRefObj, tauRefObj );
  double diff =  computeError<T>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );

  if(diff != 0.0)
  {
    printf( "orglq(): Failure Diff = %E\n", diff);
  }else{
    printf( "orglq(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}


template< typename T >
void ungqr_test()
{

  int rValC;
  int m = 128;
  int n = 256;
  int k = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aRefObj, tauRefObj;
  T *aInBuff, *aRefBuff ;
  T *tauBuff, *tauRefBuff ;
  T *work;
  char blas_uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aInBuff, aRefBuff, n*n);
  tauBuff =  new T [n-1];
  tauRefBuff =  new T [n-1];
  work = new T [lwork];

  //Call CPP function
  libflame::ungqr( &m, &n, &k, aInBuff, &m, tauBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aRefObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauRefObj );
  FLA_Obj_attach_buffer( aRefBuff, 1, m, &aRefObj );
  FLA_Obj_attach_buffer( tauRefBuff, 1, k-1, &tauRefObj ); 

  //Call C function
  rValC =  ungqr_C( aRefObj, tauRefObj );
  double diff =  computeError<T>( m, n, aRefBuff, aInBuff );
  diff +=  computeError<T>( 1, n-1, tauRefBuff, tauBuff );
  
  if(diff != 0.0)
  {
    printf( "ungqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungqr(): Success\n");
  }

  //Free up the buffers
  delete aInBuff ;
  delete tauBuff ;
  FLA_Obj_free( &aRefObj );
  FLA_Obj_free( &tauRefObj );
}
#endif

void potrf_testall_variants(){
  potrf_test<float>();
  potrf_test<double>();
  potrf_test<lapack_complex_float>();
  potrf_test<lapack_complex_double>();
}

void potf2_testall_variants(){
  potf2_test<float>();
  potf2_test<double>();
  potf2_test<lapack_complex_float>();
  potf2_test<lapack_complex_double>();
}
void getrf_testall_variants(){
  getrf_test<float>();
  getrf_test<double>();
  getrf_test<lapack_complex_float>();
  getrf_test<lapack_complex_double>();
}

void getf2_testall_variants(){
  getf2_test<float>();
  getf2_test<double>();
  getf2_test<lapack_complex_float>();
  getf2_test<lapack_complex_double>();
}

void geqrf_testall_variants(){
  geqrf_test<float>();
  geqrf_test<double>();
  geqrf_test<lapack_complex_float>();
  geqrf_test<lapack_complex_double>();
}

void geqr2_testall_variants(){
  geqr2_test<float>();
  geqr2_test<double>();
  geqr2_test<lapack_complex_float>();
  geqr2_test<lapack_complex_double>();
}

//void geqpf_testall_variants(){
//  geqpf_test<float>();
//  geqpf_test<double>();
//  geqpf_test<lapack_complex_float, float>();
//  geqpf_test<lapack_complex_double, double>();
//}

//void geqp3_testall_variants(){
//  geqp3_test<float>();
//  geqp3_test<double>();
//  geqp3_test<lapack_complex_float, float>();
//  geqp3_test<lapack_complex_double, double>();
//}


void gelqf_testall_variants(){
  gelqf_test<float>();
  gelqf_test<double>();
  gelqf_test<lapack_complex_float>();
  gelqf_test<lapack_complex_double>();
}

void gelq2_testall_variants(){
  gelq2_test<float>();
  gelq2_test<double>();
  gelq2_test<lapack_complex_float>();
  gelq2_test<lapack_complex_double>();
}

//void gelsd_testall_variants(){
//  gelsd_test<float>();
//  gelsd_test<double>();
//  gelsd_test<lapack_complex_float>();
//  gelsd_test<lapack_complex_double>();
//}

void lauum_testall_variants(){
  lauum_test<float>();
  lauum_test<double>();
  lauum_test<lapack_complex_float>();
  lauum_test<lapack_complex_double>();
}

void lauu2_testall_variants(){
  lauu2_test<float>();
  lauu2_test<double>();
  lauu2_test<lapack_complex_float>();
  lauu2_test<lapack_complex_double>();
}

void potri_testall_variants(){
  potri_test<float>();
  potri_test<double>();
  potri_test<lapack_complex_float>();
  potri_test<lapack_complex_double>();
}

void trtri_testall_variants(){
  trtri_test<float>();
  trtri_test<double>();
  trtri_test<lapack_complex_float>();
  trtri_test<lapack_complex_double>();
}

void trti2_testall_variants(){
  trti2_test<float>();
  trti2_test<double>();
  trti2_test<lapack_complex_float>();
  trti2_test<lapack_complex_double>();
}

void gehrd_testall_variants(){
  gehrd_test<float>();
  gehrd_test<double>();
  gehrd_test<lapack_complex_float>();
  gehrd_test<lapack_complex_double>();
}

void gehd2_testall_variants(){
  gehd2_test<float>();
  gehd2_test<double>();
  gehd2_test<lapack_complex_float>();
  gehd2_test<lapack_complex_double>();
}

void sytrd_testall_variants(){
  sytrd_test<float>();
  sytrd_test<double>();
}

void hetrd_testall_variants(){
  hetrd_test<lapack_complex_float, float>();
  hetrd_test<lapack_complex_double, double>();
}
  
void sytd2_testall_variants(){
  sytd2_test<float>();
  sytd2_test<double>();
}

void hetd2_testall_variants(){
  hetd2_test<lapack_complex_float, float>();
  hetd2_test<lapack_complex_double, double>();
}

void gebrd_testall_variants(){
  gebrd_test<float>();
  gebrd_test<double>();
  gebrd_test<lapack_complex_float, float>();
  gebrd_test<lapack_complex_double, double>();
}

void gebd2_testall_variants(){
  gebd2_test<float>();
  gebd2_test<double>();
  gebd2_test<lapack_complex_float, float>();
  gebd2_test<lapack_complex_double, double>();
}

void sygst_testall_variants(){
  sygst_test<float>();
  sygst_test<double>();
}

void hegst_testall_variants(){
  hegst_test<lapack_complex_float>();
  hegst_test<lapack_complex_double>();
}

void sygs2_testall_variants(){
  sygs2_test<float>();
  sygs2_test<double>();
}

void hegs2_testall_variants(){
  hegs2_test<lapack_complex_float>();
  hegs2_test<lapack_complex_double>();
}

void orgqr_testall_variants(){
  orgqr_test<float>();
  orgqr_test<double>();
  orgqr_test<lapack_complex_float>();
  orgqr_test<lapack_complex_double>();
}

void ungqr_testall_variants(){
  ungqr_test<float>();
  ungqr_test<double>();
  ungqr_test<lapack_complex_float>();
  ungqr_test<lapack_complex_double>();
}
#if 0
void orglq_testall_variants(){
  orglq_test<float>();
  orglq_test<double>();
  orglq_test<lapack_complex_float>();
  orglq_test<lapack_complex_double>();
}

void unglq_testall_variants(){
  unglq_test<float>();
  unglq_test<double>();
  unglq_test<lapack_complex_float>();
  unglq_test<lapack_complex_double>();
}

#endif
//#define Test(fnName)\
// test_ ## fnName ## (float)();
// 
 //test_ ## fnName <double>();\
 //test_ ## fnName <lapack_complex_float>();\
 //test_ ## fnName <lapack_complex_double>();

int main(int argc, char *argv[])
{
  potrf_testall_variants(); //pass
  potf2_testall_variants(); //pass
  getrf_testall_variants(); //pivot mismatch
  getf2_testall_variants();//pivot mismatch
  geqrf_testall_variants(); //pass
  geqr2_testall_variants(); //pass
  //geqpf_testall_variants(); //LAPACKE_sgeqpf not defined
  //geqp3_testall_variants(); //LAPACKE_sgeqp3 not defined
  gelqf_testall_variants(); //pass
  gelq2_testall_variants(); //pass
  //gelsd_testall_variants(); //implemenation pending
  //gelss_testall_variants(); //implemenation pending
  lauum_testall_variants(); //pass
  lauu2_testall_variants(); //pass
  potri_testall_variants(); //implemenation pending
  trtri_testall_variants(); //pass //m>128 fails
  trti2_testall_variants(); //pass //m>128 fails
  //rsyl_testall_variants(); //implemenation pending
  gehrd_testall_variants();//pass 
  gehd2_testall_variants(); //tau fails
  sytrd_testall_variants(); //pass
  hetrd_testall_variants();//pass
  sytd2_testall_variants(); //pass
  hetd2_testall_variants(); //fails //in/out buff failure--after 3 decimal point
  gebrd_testall_variants(); //pass
  gebd2_testall_variants(); //fails //alll buff mismatch-- after few decimal point
      /*Testing to be done*/
  sygst_testall_variants();//pass
  hegst_testall_variants();//pass //m>128 fails for complex float
  sygs2_testall_variants();//pass //m>64 fails for complex float
  hegs2_testall_variants();//pass //m>64 fails for complex float
  //larft
  //larfg
  //larfgp
  //orgqr_testall_variants();
  //ungqr_testall_variants();
  //ormqr
  //unmqr
  //orm2r
  //unm2r
  //orglq_testall_variants();
  //unglq_testall_variants();
}
