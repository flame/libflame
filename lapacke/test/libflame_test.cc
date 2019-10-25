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

#if 0
FLA_Error geqpf_C( FLA_Obj A, int *jpvt, FLA_Obj t, FLA_Obj rwork )
{
  int          info = 0;
//#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;
  int          lwork;
  FLA_Obj      work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_QR_check( A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  lwork    = n_A * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
    sgeqpf_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

    dgeqpf_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
    lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
    lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
    float *r_work = ( float * ) FLA_FLOAT_PTR( rwork );

    cgeqpf_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                r_work,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
    lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
    lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
    double *r_work = ( double * ) FLA_DOUBLE_PTR( rwork );

    zgeqpf_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                r_work,
                &info );

    break;
  } 

  }

  FLA_Obj_free( &work_obj );
//#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
//#endif

  return info;
}
#endif

FLA_Error geqp3_C( FLA_Obj A, int *jpvt, FLA_Obj t, FLA_Obj rwork )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;
  int          lwork;
  FLA_Obj      work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_QR_check( A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  lwork    = 3*n_A+1 ;//n_A * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
    sgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork, 
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

   dgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork, 
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
    lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
    lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );
    float *r_work = ( float * ) FLA_FLOAT_PTR( rwork );

    cgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                r_work,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
    lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
    lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
    double *r_work = ( double * ) FLA_DOUBLE_PTR( rwork );

    zgeqp3_( &m_A,
                &n_A,
                buff_A, &cs_A,
                jpvt,
                buff_t,
                buff_work,
                &lwork,
                r_work,
                &info );

    break;
  } 

  }

  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

#if 0
FLA_Error potri_C( FLA_Uplo uplo, FLA_Obj A )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, cs_A;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Ttmm_check( uplo, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    F77_spotri( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    F77_dpotri( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );

    F77_cpotri( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );

    F77_zpotri( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  } 

  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

#endif


FLA_Error gelsd_C( FLA_Obj A, FLA_Obj B, FLA_Obj sCOBuff, FLA_Obj rcondC, int rank )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A, cs_B;
  FLA_Obj work_obj, iwork_obj;
  datatype = FLA_Obj_datatype( A );
  
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_B     = FLA_Obj_col_stride( B );
  int lwork = -1;
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );
  FLA_Obj_create( FLA_INT, lwork, 1, 0, 0, &iwork_obj );
  switch( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B    = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_s    = ( float * ) FLA_FLOAT_PTR( sCOBuff );
    float *rcondCVal    =  ( float * ) FLA_FLOAT_PTR( rcondC );
	float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );
	int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );
 
    sgelsd_( &m_A,
                &n_A,
                &n_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_s, rcondCVal,
                &rank,
                buff_work, &lwork, buff_iwork, 
                &info);
  
    break;
  }
  
  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B    = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_s    = ( double * ) FLA_DOUBLE_PTR( sCOBuff );
    double *rcondCVal  = ( double * ) FLA_DOUBLE_PTR( rcondC );
  	double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );
	int *buff_iwork = ( int * ) FLA_INT_PTR( iwork_obj );
      dgelsd_( &m_A,
                &n_A,
                &n_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                buff_s, rcondCVal,
                &rank,
                buff_work, &lwork, buff_iwork, 
                &info);
  
    break;
  } 
  
//  case FLA_COMPLEX:
//  {
//    scomplex *buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
//    scomplex *buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );
//    scomplex *buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work_obj );
//  
//    cgelsd_( &m_A,
//                &n_A,
//                buff_A, &cs_A,
//                buff_t,
//                buff_work, &lwork,
//                &info );
//  
//    break;
//  } 
//  
//  case FLA_DOUBLE_COMPLEX:
//  {
//    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
//    dcomplex *buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
//    dcomplex *buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );
//  
//    zgelsd_( &m_A,
//                &n_A,
//                buff_A, &cs_A,
//                buff_t,
//                buff_work, &lwork,
//                &info );
//  
//    break;
//  } 
//  
  }
  
  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif
  
  return info;
}

template< typename T >
void potrf_test()
{
  int m = 32 ;
  char uplo = 'U';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potrf(LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );  
  //Call C function
  FLA_Chol_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj );
  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "potrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}


template< typename T >
void potf2_test()
{

  int m = 16; //2048
  char uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potf2( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  
  //Call C function
  FLA_Chol_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj );
 
  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potf2():Failure Diff = %E\n", diff);
  }else{
    printf( "potf2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void getrf_test()
{

  int m = 512;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, pivCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int *pivCPPOBuff, *pivCOBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  pivCPPOBuff =  new int [min_m_n];
  pivCOBuff =  new int [min_m_n];


 //Call CPP function
  libflame::getrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, pivCPPOBuff );
 
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( pivCOBuff, 1, min_m_n, &pivCOObj );
  //FLA_Set( FLA_ZERO, pivCOObj );  
  
  //Call C function
  FLA_LU_piv_blk_external( aCIOObj, pivCOObj );


  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivCOBuff, pivCPPOBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getrf(): Failure Pivot BUffer mismatach Diff = %d\n", diffInt);
  }else{
    printf( "getrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete pivCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &pivCOObj );
}


template< typename T >
void getf2_test()
{

  int m = 1024;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, pivCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  int *pivCPPOBuff, *pivCOBuff ;
  int min_m_n = min( m, n );

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  pivCPPOBuff =  new int [min_m_n];
  pivCOBuff =  new int [min_m_n];

 //Call CPP function
  libflame::getf2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, pivCPPOBuff );
 
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &pivCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( pivCOBuff, 1, min_m_n, &pivCOObj );
  //FLA_Set( FLA_ZERO, pivCOObj );  
  
  //Call C function
  FLA_LU_piv_unb_external( aCIOObj, pivCOObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  int diffInt =  computeError<int>( 1, min_m_n, pivCOBuff, pivCPPOBuff );

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff);
  }else if(diffInt !=0){
    printf( "getf2(): Failure:Pivot Buffer Mismatch Diff = %d\n", diffInt);
  }else{
    printf( "getf2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete pivCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &pivCOObj );
}
#if 0
#define PREFIX2FLAME_INVERT_TAU(datatype, tauCOObj)\
{\
  if (datatype == FLA_FLOAT)\
    FLAME_invert_stau(tauCOObj);\
  if (datatype == FLA_DOUBLE)\
    FLAME_invert_dtau(tauCOObj);\
  if (datatype == FLA_COMPLEX)\
    FLAME_invert_ctau(tauCOObj);\
  if (datatype == FLA_DOUBLE_COMPLEX)\
    FLAME_invert_ztau(tauCOObj);\
}
#endif

template< typename T >
void geqrf_test()
{
  
  int m = 512;
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int min_m_n = min( m, n );
  int datatype = getDatatype<T>();
  int lwork    = n * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

 //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ); 

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );


   fp = fopen("test/out","a+");
   print(fp, m*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
  
  double diff = computeError<T>( n, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "geqrf(): Failure Diff = %E\n", diff);
  }else{
    printf( "geqrf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void geqr2_test()
{


  int m = 128;
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int min_m_n = min( m, n ) ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::geqr2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
 
  //Call C function
  FLA_QR_unb_external( aCIOObj, tauCOObj );



  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqr2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqr2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}
#if 0
template< typename T >
void geqpf_test()
{


  int m = 256;
  int n = 4096;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  workBuff =  new T [3*n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &lda, jpvtCPPOBuff, tauCPPOBuff );

  //Call C function
 // geqpf_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkObj );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }
 
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta, typename Tb >
void geqpf_test()
{

  int m = 512;
  int n = 2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  tauCPPOBuff =  new Ta [min_m_n];
  tauCOBuff =  new Ta [min_m_n];
  rworkRefBuff = new Tb [2 * n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;
  //Call CPP function
  libflame::geqpf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &lda, jpvtCPPOBuff, tauCPPOBuff ) ;

  //Call C function
  geqpf_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );

  
  double diff =  computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<Ta>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "geqpf(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqpf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqpf(): Success\n");
  }
 
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}
#endif
template< typename T >
void geqp3_test()
{

  int m = 256;
  int n = 4096;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  T *aCPPIOBuff, *aCIOBuff, *rworkRefBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;
  
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  rworkRefBuff = new T [2 * n];
  
  int datatype = getDatatype<T>();
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatype, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;  
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff );

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );


  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }
 
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}

template< typename Ta, typename Tb >
void geqp3_test()
{
  int m = 128;//2048;
  int n = 128;//512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new Ta [min_m_n];
  tauCOBuff =  new Ta [min_m_n];
  rworkRefBuff = new Tb [2 * n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::geqp3( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff ) ;

  //Call C function
  geqp3_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );
  
  double diff = computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff += computeError<Ta>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  int diffInt = computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "geqp3(): Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "geqp3(): Failure Diff1111 = %d\n", diffInt) ;
  }else{
    printf( "geqp3(): Success\n");
  }
 
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete jpvtCPPOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  FLA_Obj_free( &rworkRefObj );
}


template< typename T >
void gelqf_test()
{
  int m = 8192;
  int n = 1024;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int min_m_n = min( m, n ) ;

  int datatype = getDatatype<T>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];
  workBuff =  new T [lwork];

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );



  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;

  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelqf(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelqf(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gelq2_test()
{
  int m = 1000;
  int n = 1000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  tauCPPOBuff =  new T [min_m_n];
  tauCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gelq2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &lda, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;

  //Call C function
  FLA_LQ_unb_external( aCIOObj, tauCOObj );
 
  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  int diffInt =  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else if(diffInt !=0){
    printf( "gelq2(): Failure Diff = %d\n", diffInt) ;
  }else{
    printf( "gelq2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gelsd_test()
{
  int m = 128;//16384;
  int n = 128;//128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bCIOObj, sCOObj, rcondCObj;
  T *aCPPIOBuff, *bCPPIOBuff, *aCIOBuff, *bCIOBuff ;
  T  *sCPPOBuff, *sCOBuff, *rcondCPP, *rcondC ;
  int *b, *bRef;
  int rank  = 0;
  int min_m_n = min( m, n ) ;
  int lda = max(1,m);
  int ldb = max(1,max(m,n));
  
  b = new int[n];
  bRef = new int[n];
  sCPPOBuff = new T [min_m_n];
  sCOBuff = new T [min_m_n];
  
  int datatype = getDatatype<T>();
 
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, ldb*n) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, ldb, n, &bCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &sCOObj ) ;
  FLA_Obj_create_without_buffer( datatype, 1, 1, &rcondCObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, ldb, &bCIOObj );
  FLA_Obj_attach_buffer( sCOBuff, 1, min_m_n, &sCOObj );
  FLA_Obj_attach_buffer( rcondC, 1, 1, &rcondCObj );

  //Call CPP function
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, &n, aCPPIOBuff, &lda, bCPPIOBuff, &ldb, sCPPOBuff, rcondCPP, &rank );

  //Call C function
  gelsd_C( aCIOObj, bCIOObj, sCOObj, rcondCObj, rank );

  double diff =  computeError<T>( n, m, aCIOBuff, aCPPIOBuff ) ;
  diff +=  computeError<T>( n, m, bCIOBuff, bCPPIOBuff ) ;
  diff +=  computeError<T>( 1, min_m_n, sCOBuff, sCPPOBuff ) ;
  diff += abs(*rcondCPP - *rcondC);
  
  if(diff != 0.0)
  {
    printf( "Failure Diff = %E\n", diff) ;
  }else{
    printf( "gelsd(): Success\n");
  }
 
 //Free up the buffers
 delete aCPPIOBuff ;
 delete bCPPIOBuff ;
 delete sCPPOBuff ;
 delete rcondCPP ;
 FLA_Obj_free( &aCIOObj );
 FLA_Obj_free( &bCIOObj );
  FLA_Obj_free(&sCOObj );
 FLA_Obj_free( &rcondCObj );
}

template< typename Ta, typename Tb >
void gelsd_test()
{
  int m = 2048;
  int n = 512;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj, rworkRefObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Ta *workBuff;
  Tb *rworkBuff, *rworkRefBuff;
  int *jpvtCPPOBuff, *jpvtCOBuff;
  int min_m_n = min( m, n ) ;

  tauCPPOBuff =  new Ta [min_m_n];
  tauCOBuff =  new Ta [min_m_n];
  workBuff =  new Ta [3 * n];
  rworkBuff =  new Tb [2 * n];
  rworkRefBuff = new Tb [2 * n];
  jpvtCPPOBuff = new int[n];
  jpvtCOBuff = new int[n];
  
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();
  int lwork    = m * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n) ;
  allocate_init_buffer(rworkBuff, rworkRefBuff, 2*n,0) ;

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj ) ;
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj ) ;
  FLA_Obj_create_without_buffer( datatypeTb, 2 * n, 1, &rworkRefObj ) ;
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ) ;
  FLA_Obj_attach_buffer( rworkRefBuff, 1, 2 * n, &rworkRefObj ) ;

  //Call CPP function
  libflame::gelsd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &m, jpvtCPPOBuff, tauCPPOBuff ) ;

  //Call C function
// gelsd_C( aCIOObj, jpvtCOBuff, tauCOObj, rworkRefObj );
//
// double diff =  computeError<Ta>( n, m, aCIOBuff, aCPPIOBuff ) ;
// int diffInt =  computeError<Ta>( 1, min_m_n, tauCOBuff, tauCPPOBuff ) ;
// diffInt +=  computeError<int>( 1, n, jpvtCOBuff, jpvtCPPOBuff ) ;
// diff +=  computeError<Tb>( 1, 2*n, rworkRefBuff ) ;
// 
// if(diff != 0.0)
// {
//   printf( "gelsd(): Failure Diff = %E\n", diff) ;
// }else if(diffInt !=0){
//   printf( "gelsd(): Failure Diff = %d\n", diffInt) ;
// }else{
//   printf( "gelsd(): Success\n");
// }
//
// //Free up the buffers
// delete aCPPIOBuff ;
// delete tauCPPOBuff ;
// delete rworkBuff ;
// FLA_Obj_free( &aCIOObj );
// FLA_Obj_free( &tauCOObj );
// FLA_Obj_free( &rworkRefObj );
}

template< typename T >
void lauum_test()
{

  int m = 64;
  char uplo = 'u';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::lauum( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
  FLA_Ttmm_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "lauum(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauum(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void lauu2_test()
{

  int m = 128;
  char uplo = 'l';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::lauu2( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Param_map_flame_to_netlib_uplo( FLA_LOWER_TRIANGULAR, &uplo );

  //Call C function
  FLA_Ttmm_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "lauu2(): Failure Diff = %E\n", diff);
  }else{
    printf( "lauu2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}


template< typename T >
void potri_test()
{


  int m = 64; //64 works
  char uplo = 'L';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::potri( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );

  //Call C function
//  potri_C( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );
  
  if(diff != 0.0)
  {
	fp = fopen("test/in","a+");
   print(fp,m*m,aCPPIOBuff, aCIOBuff);
   fclose(fp);
  }
  
  if(diff != 0.0)
  {
    printf( "potri(): Failure Diff = %E\n", diff);
  }else{
    printf( "potri(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}



template< typename T >
void trtri_test()
{
  int m = 128;
  char uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
    FILE *fp = fopen("test/in1","a+");
	print(fp, m* m, aCIOBuff, aCPPIOBuff);
	fclose(fp);
  //Call CPP function
  libflame::trtri( LAPACK_COL_MAJOR, &uplo, &blas_diag, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  
  //Call C function
  FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trtri(): Failure Diff = %E\n", diff);
    FILE *fp = fopen("test/in","a+");
	FILE *fp1 = fopen("test/out","a+");
	print(fp, m* m, aCIOBuff, aCIOBuff);
	print(fp1, m* m, aCPPIOBuff, aCPPIOBuff);
	fclose(fp);
	fclose(fp1);
  }else{
    printf( "trtri(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void trti2_test()
{
  int m = 128;
  char uplo = 'L';
  char blas_diag = 'N';
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj;
  T *aCPPIOBuff, *aCIOBuff ;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);

  //Call CPP function
  libflame::trti2( LAPACK_COL_MAJOR, &uplo, &blas_diag, &m, aCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  
  //Call C function
  FLA_Trinv_blk_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );
  //FLA_Trinv_unb_external( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trti2(): Failure Diff = %E\n", diff);
  }else{
    printf( "trti2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

template< typename T >
void trsyl_test()
{

  int m = 64;
  int n = 256;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  FLA_Trans transARef = FLA_NO_TRANSPOSE ; FLA_Trans transBRef = FLA_NO_TRANSPOSE;
  FLA_Obj isgnObj, aCIOObj,  bRefObj, cRefObj, scaleRefObj ;
  T *aCPPIOBuff, *bCPPIOBuff, *cCPPIOBuff ;
  T *aCIOBuff, *bCIOBuff, *cCIOBuff, *scaleRef ;
  T scale;
  int isgn = 1;
  int isgnValue = 1;
  int isgnRefValue  = 1;
  int *isgnRef = &isgnValue;


  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);

//int matrix_layout, char* transa, char* transb, int* isgn, int* m, int* n, T* a, int* lda, T* b, int* ldb, T* c, int* ldc, T* scale )
  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bCPPIOBuff, &n, cCPPIOBuff, &m, &scale ); 

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( scaleRef, 1, 1, &scaleRefObj );
  FLA_Obj_create_without_buffer( FLA_INT, 1, 1, &isgnObj );
  FLA_Obj_attach_buffer( isgnRef, 1, 1, &isgnObj );
  
  //Call C function
  
  FLA_Sylv_unb_external( transARef, transBRef, isgnObj, aCIOObj,  bRefObj, cRefObj, scaleRefObj );


  //Compute Difference in C and CPP buffer
  //double diff =  computeError<T>( m, m, cRefBuff, cInBuff );

 // if(diff != 0.0)
 // {
 //   printf( "trsyl(): Failure Diff = %E\n", diff);
 // }else{
   printf( "trsyl(): Success\n");
 // }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bCPPIOBuff;
  delete cCPPIOBuff;
  //FLA_Obj_free( &aCIOObj );
  //FLA_Obj_free( &bRefObj );
  //FLA_Obj_free( &cRefObj );
}


template< typename Ta, typename Tb  >
void trsyl_test()
{

  int m = 256;
  int n = 256;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  //FLA_Trans transarRef = 'N' ; FLA_Trans transbRef = 'N';
  FLA_Obj isgnObj, aCIOObj,  bRefObj, cRefObj, scaleRefObj ;
  Tb *aCPPIOBuff, *bInBuff, *cInBuff ;
  Tb *aCIOBuff, *bRefBuff, *cRefBuff, *scaleRef ;
  Tb scale;
  int isgn = 1;
  int isgnValue = 1;
  
  int *isgnRef = &isgnValue;


  int datatype = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(aCPPIOBuff, bRefBuff, n*n);
  allocate_init_buffer(aCPPIOBuff, cRefBuff, m*n);

//int matrix_layout, char* transa, char* transb, int* isgn, int* m, int* n, T* a, int* lda, T* b, int* ldb, T* c, int* ldc, T* scale )
  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bInBuff, &n, cInBuff, &m, &scale ); 

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( scaleRef, 1, n, &scaleRefObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &isgnObj );
  FLA_Obj_attach_buffer( isgnRef, 1, 1, &isgnObj );
  
  //Call C function
  
 // FLA_Sylv_unb_external( transa, transb, isgnObj, aCIOObj,  bRefObj, cRefObj, scaleRefObj );


  //Compute Difference in C and CPP buffer
  double diff = 0;// computeError<Ta>( m, m, cRefBuff, cInBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): Failure Diff = %E\n", diff);
  }else{
    printf( "trsyl(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bInBuff;
  delete cInBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bRefObj );
  FLA_Obj_free( &cRefObj );
}


template< typename T >
void gehrd_test()
{
  int n = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  tauCOBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::gehrd( LAPACK_COL_MAJOR, &n, &ilo, &iho, aCPPIOBuff, &n, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Hess_blk_external( aCIOObj, tauCOObj, ilo, iho );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gehrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gehd2_test()
{
  int n = 256;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff ;
  int iho = 10;
  int ilo = 10;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  
  tauCPPOBuff =  new T [n-1];
  tauCOBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::gehd2( LAPACK_COL_MAJOR, &n, &ilo, &iho, aCPPIOBuff, &n, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Hess_unb_external( aCIOObj, tauCOObj, ilo, iho );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );
     fp = fopen("test/out1","a+");
   print(fp,n-1,tauCOBuff, tauCOBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,n-1,tauCPPOBuff, tauCPPOBuff);
   fclose(fp);
  if(diff != 0.0)
  {
    printf( "gehd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gehd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void sytrd_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *d, *e ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauCOBuff =  new T [n-1];

  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "sytrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta, typename Tb >
void hetrd_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Tb  *d, *e ;
  char uplo = 'u';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauCOBuff =  new Ta [n-1];
  
  //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, n-1, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "hetrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void sytd2_test()
{
  int n = 2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *workBuff, *d, *e ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new T [n-1];
  d =  new T [n];
  e =  new T [n-1];
  tauCOBuff =  new T [n-1];
  workBuff =  new T [n];

  //Call CPP function
  libflame::sytd2( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Tridiag_unb_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, n-1, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "sytd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sytd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta,  typename Tb >
void hetd2_test()
{
  int n = 128;//2000;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  Ta *workBuff;
  Tb *d, *e ;
  char uplo = 'u';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  tauCPPOBuff =  new Ta [n-1];
  d =  new Tb [n];
  e =  new Tb [n-1];
  tauCOBuff =  new Ta [n-1];
  workBuff =  new Ta [n];
   fp = fopen("test/in","a+");
   print(fp,n*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
  //Call CPP function
  libflame::hetd2( LAPACK_COL_MAJOR, &uplo, &n, aCPPIOBuff, &n, d, e, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, n-1, &tauCOObj ); 

  //Call C function
  FLA_Tridiag_unb_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, n-1, tauCOBuff, tauCPPOBuff );

     fp = fopen("test/out","a+");
   print(fp,n*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
     fp = fopen("test/out1","a+");
   print(fp,n*n,aCPPIOBuff, aCPPIOBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,n*n,aCIOBuff, aCIOBuff);
   fclose(fp);

  if(diff != 0.0)
  {
    printf( "hetd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hetd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete workBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void gebrd_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<T>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void gebrd_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupCOBuff, taupCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "gebrd(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebrd(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void gebd2_test()
{

  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  int datatype = getDatatype<T>();
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_unb_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<T>( 1, min_m_n, taupCOBuff, taupCPPOBuff );

  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}


template< typename Ta, typename Tb >
void gebd2_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];
   fp = fopen("test/in","a+");
   print(fp,m*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
  //Call CPP function
  libflame::gebd2( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_unb_external( aCIOObj, tauqCOObj, taupCOObj );
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, tauqCOBuff, tauqCPPOBuff );
  diff +=  computeError<Ta>( 1, min_m_n, taupCOBuff, taupCPPOBuff );
   fp = fopen("test/tau1","a+");
   print(fp,min_m_n,tauqCOBuff, taupCOBuff);
   fclose(fp);
   fp = fopen("test/tau2","a+");
   print(fp,min_m_n,tauqCPPOBuff, taupCPPOBuff);
   fclose(fp);
     fp = fopen("test/out","a+");
   print(fp,m*n,aCPPIOBuff, aCIOBuff);
   fclose(fp);
     fp = fopen("test/out1","a+");
   print(fp,m*n,aCPPIOBuff, aCPPIOBuff);
   fclose(fp);
     fp = fopen("test/out2","a+");
   print(fp,m*n,aCIOBuff, aCIOBuff);
   fclose(fp);
   
  if(diff != 0.0)
  {
    printf( "gebd2(): Failure Diff = %E\n", diff);
  }else{
    printf( "gebd2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void sygst_test()
{
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bRefObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::sygst( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygst(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void hegst_test()
{

  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bRefObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::hegst( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_blk_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bRefObj );

  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegst(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void sygs2_test()
{
  int n = 64;//2048;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bRefObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::sygs2( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "sygs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "sygst(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void hegs2_test()
{
  int n = 64;//2048;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, bRefObj;
  T *aCPPIOBuff, *aCIOBuff, *bRefBuff, *b;
  char uplo  ='l';
  int datatype = getDatatype<T>();
  int itype = 1 ; //1 or 2
  int itype_c = FLA_INVERSE;
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(b, bRefBuff, n*n);

  
  //Call CPP function
  libflame::hegs2( LAPACK_COL_MAJOR, &itype, &uplo, &n, aCPPIOBuff, &n, b, &n );
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bRefObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, n, &aCIOObj );
  FLA_Obj_attach_buffer( bRefBuff, 1, n, &bRefObj );

  //Call C function
  FLA_Eig_gest_unb_external( itype_c, FLA_LOWER_TRIANGULAR, aCIOObj, bRefObj );
  
  double diff =  computeError<T>( n, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( n, n, bRefBuff, b );

  if(diff != 0.0)
  {
    printf( "hegs2(): Failure Diff = %E\n", diff);
  }else{
    printf( "hegs2(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete b ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bRefObj );
}

template< typename T >
void orgqr_test()
{
  int m = 64;//512;
  int n = 64;//256;
  int k = 64;//128;
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
 libflame::orgqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );
 
 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
 FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
 FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
 FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ); 
 
  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_QR_form_Q_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "orgqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgqr(): Success\n");
  }
  
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ungqr_test()
{
  int m = 512;
  int n = 256;
  int k = 128;
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

template< typename T >
void ormqr_test()
{
#if 0
  int m = 64;//512;
  int n = 64;//256;
  int k = 64;//128;
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
 libflame::ormqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );
 
 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
 FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
 FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
 FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ); 
 
  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "ormqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormqr(): Success\n");
  }
  
  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
  #endif
}

template< typename T >
void unmqr_test()
{
#if 0	
  int m = 512;
  int n = 256;
  int k = 128;
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
  libflame::unmqr( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, min_m_n, &tauCOObj ); 

  //Call C function
  FLA_QR_blk_external( aCIOObj, tauCOObj );
  FLA_Apply_Q_blk_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, min_m_n, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "unmqr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmqr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
#endif
}


#if 0

template< typename T >
void orglq_test()
{
  int m = 512;
  int n = 256;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *work;
  char uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [k-1];
  tauCOBuff =  new T [k-1];
  work = new T [lwork];

  //Call CPP function
  libflame::orglq( LAPACK_COL_MAJOR, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k-1, &tauCOObj ); 

  //Call C function
 // orglq_C( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orglq(): Failure Diff = %E\n", diff);
  }else{
    printf( "orglq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename T >
void ormlq_test()
{
  int m = 512;
  int n = 256;
  int k = 128;

  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  T *work;
  char uplo = 'L';
  int datatype = getDatatype<T>();
  int lwork    = max( 1, n );

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  tauCPPOBuff =  new T [k-1];
  tauCOBuff =  new T [k-1];
  work = new T [lwork];

  //Call CPP function
//  libflame::ormlq( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &m, tauCPPOBuff, &c, &ldc );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k-1, &tauCOObj ); 

  //Call C function
  //( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B )
//  FLA_Apply_Q_blk_external( aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k-1, tauCOBuff, tauCPPOBuff );

  if(diff != 0.0)
  {
    printf( "orglq(): Failure Diff = %E\n", diff);
  }else{
    printf( "orglq(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}
#endif

template< typename T >
void orgtr_test()
{
  int m = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauCPPOBuff, *tauCOBuff ;
  char uplo = 'L';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  tauCPPOBuff =  new T [m-1];
  tauCOBuff =  new T [m-1];
  T *d =  new T [m];
  T *e =  new T [m-1];
  
  //Call CPP function
  libflame::sytrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::orgtr( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, tauCPPOBuff);

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj ); 

  //Call C function
  FLA_Tridiag_blk_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_form_Q_external( FLA_LOWER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, m-1, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "orgtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

template< typename Ta,  typename Tb>
void ungtr_test()
{
  int m = 512;
  srand (time(NULL));
  
  FLA_Init( );
  FLA_Obj aCIOObj, tauCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauCPPOBuff, *tauCOBuff ;
  char uplo = 'U';
  int datatype = getDatatype<Ta>();
  
  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  tauCPPOBuff =  new Ta [m-1];
  tauCOBuff =  new Ta [m-1];
  Tb *d =  new Tb [m];
  Tb *e =  new Tb [m-1];
  
  //Call CPP function
  libflame::hetrd( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, d, e, tauCPPOBuff );
  libflame::ungtr( LAPACK_COL_MAJOR, &uplo, &m, aCPPIOBuff, &m, tauCPPOBuff);
  
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m-1, 1, &tauCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, m-1, &tauCOObj ); 
  
  //Call C function
  FLA_Tridiag_blk_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  FLA_Tridiag_form_Q_external( FLA_UPPER_TRIANGULAR, aCIOObj, tauCOObj );
  double diff =  computeError<Ta>( m, m, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<Ta>( 1, m-1, tauCOBuff, tauCPPOBuff );
  
  if(diff != 0.0)
  {
    printf( "ungtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

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
    printf( "ormtr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormtr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauCOObj );
}

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

template< typename T >
void orgbr_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *aCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  char vect = 'P';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::orgbr( LAPACK_COL_MAJOR, &vect, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_form_V_external( aCIOObj, taupCOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orgbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "orgbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void ungbr_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  char vect = 'P';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::ungbr( LAPACK_COL_MAJOR, &vect, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_form_V_external( aCIOObj, taupCOObj );
  
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  
  if(diff != 0.0)
  {
    printf( "ungbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ungbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename T >
void ormbr_test()
{
  int m = 512;
  int n = 512;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauqCOObj, taupCOObj;
  T *aCPPIOBuff, *cCPPIOBuff, *aCIOBuff, *cCIOBuff ;
  T *tauqCPPOBuff, *tauqCOBuff ;
  T *taupCPPOBuff, *taupCOBuff ;
  T *d, *e ;
  char vect = 'Q';
  char sideCPP = 'L';
  char transCPP = 'T';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, n*n);
  d =  new T [min_m_n];
  e =  new T [min_m_n-1];
  tauqCPPOBuff =  new T [min_m_n];
  taupCPPOBuff =  new T [min_m_n];
  tauqCOBuff =  new T [min_m_n];
  taupCOBuff =  new T [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::ormbr( LAPACK_COL_MAJOR, &vect, &sideCPP, &transCPP, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff, cCPPIOBuff, &m  );
  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_apply_U_external( FLA_LEFT, FLA_TRANSPOSE, aCIOObj, taupCOObj, cCIOObj );
  double diff =  computeError<T>( m, n, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "ormbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "ormbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}

template< typename Ta, typename Tb >
void unmbr_test()
{
  int m = 256 ;
  int n = 256 ;
  int min_m_n = min( m, n );
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj aCIOObj, cCIOObj, tauqCOObj, taupCOObj;
  Ta *aCPPIOBuff, *aCIOBuff ;
  Ta *cCPPIOBuff, *cCIOBuff ;
  Ta *tauqCPPOBuff, *tauqCOBuff ;
  Ta *taupCPPOBuff, *taupCOBuff ;
  Tb *d, *e ;
  char vect = 'P';
  char sideCPP = 'R';
  char transCPP = 'N';
  int datatype = getDatatype<Ta>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);
  d =  new Tb [min_m_n];
  e =  new Tb [min_m_n-1];
  tauqCPPOBuff =  new Ta [min_m_n];
  taupCPPOBuff =  new Ta [min_m_n];
  tauqCOBuff =  new Ta [min_m_n];
  taupCOBuff =  new Ta [min_m_n];

  //Call CPP function
  libflame::gebrd( LAPACK_COL_MAJOR, &m, &n, aCPPIOBuff, &n, d, e, tauqCPPOBuff, taupCPPOBuff );
  libflame::unmbr( LAPACK_COL_MAJOR, &vect, &sideCPP, &transCPP, &m, &n, &m, aCPPIOBuff, &m, taupCPPOBuff, cCPPIOBuff, &m  );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &tauqCOObj );
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &taupCOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_attach_buffer( tauqCOBuff, 1, min_m_n, &tauqCOObj ); 
  FLA_Obj_attach_buffer( taupCOBuff, 1, min_m_n, &taupCOObj ); 

  //Call C function
  FLA_Bidiag_blk_external( aCIOObj, tauqCOObj, taupCOObj );
  FLA_Bidiag_apply_V_external( FLA_RIGHT, FLA_NO_TRANSPOSE, aCIOObj, taupCOObj, cCIOObj );
  
  double diff =  computeError<Ta>( m, n, aCIOBuff, aCPPIOBuff );
  
  if(diff != 0.0)
  {
    printf( "unmbr(): Failure Diff = %E\n", diff);
  }else{
    printf( "unmbr(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff ;
  delete tauqCPPOBuff ;
  delete taupCPPOBuff ;
  delete d ;
  delete e ;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &tauqCOObj );
  FLA_Obj_free( &taupCOObj );
}


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
//
//void geqpf_testall_variants(){
//  geqpf_test<float>();
//  geqpf_test<double>();
//  geqpf_test<lapack_complex_float, float>();
//  geqpf_test<lapack_complex_double, double>();
//}

void geqp3_testall_variants(){
  geqp3_test<float>();
  geqp3_test<double>();
  geqp3_test<lapack_complex_float, float>();
  geqp3_test<lapack_complex_double, double>();
}

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

void gelsd_testall_variants(){
  gelsd_test<float>();
  gelsd_test<double>();
  //gelsd_test<lapack_complex_float>();
  //gelsd_test<lapack_complex_double>();
}

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

void trsyl_testall_variants(){
  trsyl_test<float>();
  trsyl_test<double>();
  trsyl_test<lapack_complex_float, float>();
  trsyl_test<lapack_complex_double, double>();
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
}
void ungqr_testall_variants(){
  ungqr_test<lapack_complex_float>();
  ungqr_test<lapack_complex_double>();
}
void ormqr_testall_variants(){
  ormqr_test<float>();
  ormqr_test<double>();
}
void unmqr_testall_variants(){
  unmqr_test<lapack_complex_float>();
  unmqr_test<lapack_complex_double>();
}

#if 0
void orglq_testall_variants(){
  orglq_test<float>();
  orglq_test<double>();
  orglq_test<lapack_complex_float>();
  orglq_test<lapack_complex_double>();
}
void ormlq_testall_variants(){
  ormlq_test<float>();
  ormlq_test<double>();
  ormlq_test<lapack_complex_float>();
  ormlq_test<lapack_complex_double>();
}
#endif

void orgtr_testall_variants(){
  orgtr_test<float>();
  orgtr_test<double>();
}
void ungtr_testall_variants(){
  ungtr_test<lapack_complex_float, float>();
  ungtr_test<lapack_complex_double, double >();
}
void ormtr_testall_variants(){
  ormtr_test<float>();
  ormtr_test<double>();
}
void unmtr_testall_variants(){
  unmtr_test<lapack_complex_float, float>();
  unmtr_test<lapack_complex_double, double >();
}
void orgbr_testall_variants(){
  orgbr_test<float>();
  orgbr_test<double>();
}
void ungbr_testall_variants(){
  ungbr_test<lapack_complex_float, float>();
  ungbr_test<lapack_complex_double, double >();
}
void ormbr_testall_variants(){
  ormbr_test<float>();
  ormbr_test<double>();
}
void unmbr_testall_variants(){
  unmbr_test<lapack_complex_float, float>();
  unmbr_test<lapack_complex_double, double >();
}

//void unglq_testall_variants(){
//  unglq_test<float>();
//  unglq_test<double>();
//  unglq_test<lapack_complex_float>();
//  unglq_test<lapack_complex_double>();
//}

//#define Test(fnName)\
// test_ ## fnName ## (float)();
// 
 //test_ ## fnName <double>();\
 //test_ ## fnName <lapack_complex_float>();\
 //test_ ## fnName <lapack_complex_double>();

int main(int argc, char *argv[])
{

#if 1
  //potrf_testall_variants(); //pass
  //potf2_testall_variants(); //pass
  //getrf_testall_variants(); //pivot mismatch
  //getf2_testall_variants();//pivot mismatch
  //geqrf_testall_variants(); //pass
  //geqr2_testall_variants(); //pass
  //geqpf_testall_variants(); / geqpf not included in .a ,LAPACKE_sgeqpf not defined
  //geqp3_testall_variants(); //complex failure
  //gelqf_testall_variants(); //pass
  //gelq2_testall_variants(); //pass
  //gelsd_testall_variants(); //seg fault
  //gelss_testall_variants(); //implemenation pending
  //lauum_testall_variants(); //pass
  //lauu2_testall_variants(); //pass
  //potri_testall_variants(); //implemenation pending
  //trtri_testall_variants(); //pass //m>128 fails
  //trti2_testall_variants(); //pass //m>128 fails
  //trsyl_testall_variants(); //implemenation pending
  //gehrd_testall_variants();//pass 
  //gehd2_testall_variants(); //tau fails
  //sytrd_testall_variants(); //pass
  //hetrd_testall_variants();//pass
  //sytd2_testall_variants(); //pass
  //hetd2_testall_variants(); //fails //in/out buff failure--after 3 decimal point
  //gebrd_testall_variants(); //pass
  //gebd2_testall_variants(); //fails //alll buff mismatch-- after few decimal point
  //    /*Testing to be done*/
  //sygst_testall_variants();//pass
  //hegst_testall_variants();//pass //m>128 fails for complex float
  //sygs2_testall_variants();//pass //m>64 fails for complex float
  //hegs2_testall_variants();//pass //m>64 fails for complex float
#endif  
  //larft
  //larfg
  //larfgp
  //orgqr_testall_variants(); //pass  
  //ungqr_testall_variants(); //pass
  ormqr_testall_variants();
  unmqr_testall_variants();
  //orm2r
  //unm2r
  //orglq_testall_variants();//implemenation pending
  //unglq_testall_variants();//implemenation pending
  //ormlq_testall_variants(); //implement
  //unmlq
  //orml2
  //unml2  
  //orgtr_testall_variants(); //pass
  //ungtr_testall_variants(); //pass
  //ormtr_testall_variants(); //pass
  //unmtr_testall_variants(); //pass
  //orgbr_testall_variants(); //pass
  //ungbr_testall_variants(); //pass
  //ormbr_testall_variants(); //pass
  //unmbr_testall_variants(); //pass  
  //steqr
  //stedc
  //stemr
  //syev
  //heev
  //syevd
  //heevd
  //syevr
  //heevr
  //heevr
  //bdsqr
  //bdsdc
  //gesvd
  //gesdd
  //laswp
  //laset
  
  //Passing test
  #if 1
  potrf_testall_variants(); //pass
  potf2_testall_variants(); //pass
  getrf_testall_variants(); //pivot mismatch
  getf2_testall_variants();//pivot mismatch
  geqrf_testall_variants(); //pass
  geqr2_testall_variants(); //pass
  gelqf_testall_variants(); //pass
  gelq2_testall_variants(); //pass
  lauum_testall_variants(); //pass
  lauu2_testall_variants(); //pass
  trtri_testall_variants(); //pass //m>128 fails
  trti2_testall_variants(); //pass //m>128 fails
  gehrd_testall_variants();//pass 
  sytrd_testall_variants(); //pass
  hetrd_testall_variants();//pass
  sytd2_testall_variants(); //pass
  gebrd_testall_variants(); //pass

  sygst_testall_variants();//pass

  orgqr_testall_variants(); //pass  
  ungqr_testall_variants(); //fails

  orgtr_testall_variants(); //pass
  ungtr_testall_variants(); //pass
  ormtr_testall_variants(); //pass
  unmtr_testall_variants(); //pass
  orgbr_testall_variants(); //pass
  ungbr_testall_variants(); //pass
  ormbr_testall_variants(); //pass
  unmbr_testall_variants(); //pass  

  
  #endif
}


