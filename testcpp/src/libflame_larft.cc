
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

/*! @file liblame_larft.cc
 *  libflame_larft.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

FLA_Error larft_C(FLA_Direct direct, FLA_Direct storev, FLA_Obj A, FLA_Obj tauObj, FLA_Obj T)
{
  int          info = 0;

#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  datatype = FLA_Obj_datatype( A );

  int m_A      = FLA_Obj_length( A );
  int n_A      = FLA_Obj_width( A );
  int cs_A     = FLA_Obj_col_stride( A );
  int ldt      = FLA_Obj_length( T );
  char char_direct,  char_storev;
  FLA_Param_map_flame_to_netlib_direct( direct, &char_direct );
  FLA_Param_map_flame_to_netlib_storev( storev, &char_storev );
  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_T    = ( float * ) FLA_FLOAT_PTR( T );
      float *buff_tau    = ( float * ) FLA_FLOAT_PTR( tauObj );
      slarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_T    = ( double * ) FLA_DOUBLE_PTR( T );
      double *buff_tau    = ( double * ) FLA_DOUBLE_PTR( tauObj );
      dlarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_T    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( T );
      lapack_complex_float *buff_tau    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( tauObj );
      clarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_T    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( T );
      lapack_complex_double *buff_tau    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( tauObj );
      zlarft_( &char_direct, &char_storev, &m_A, &n_A, buff_A, &cs_A, buff_tau, buff_T, &ldt);
      break;
    }
  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

template< typename T >
void larft_test()
{
  int n = 64;//256;
  int k = 64;//128;
 srand (time(NULL));

 FLA_Init( );
 FLA_Obj vCIObj, tauCIObj, tCOObj;
 T *vCPPIBuff, *vCIBuff, *tCPPOBuff, *tCOBuff ;
 T *tauCPPIBuff, *tauCIBuff ;
 int datatype = getDatatype<T>();
 char direct  = 'F';
 char storev = 'C';
 int ldv = n;
 int ldt = k;
 //Allocate and initialize buffers for C and CPP functions with random values
 allocate_init_buffer(vCPPIBuff, vCIBuff, ldv*k);
 tauCPPIBuff =  new T [k];
 tauCIBuff =  new T [k];
 tCPPOBuff =  new T [ldt*k];
 tCOBuff =  new T [ldt*k];
 for(int i =0; i <ldt*k; i++)
 {
    tCPPOBuff[i] = 0;
    tCOBuff[i] = 0;
 }

 //Call CPP function
  libflame::geqrf( LAPACK_COL_MAJOR, &n, &k, vCPPIBuff, &n, tauCPPIBuff );
  libflame::larft( LAPACK_COL_MAJOR, &direct, &storev, &n, &k, vCPPIBuff, &ldv, tauCPPIBuff, tCPPOBuff, &ldt );

 //Allocate Object for C function and copy already allocated and filled buffer
 FLA_Obj_create_without_buffer( datatype, ldv, k, &vCIObj );
 FLA_Obj_create_without_buffer( datatype, k, 1, &tauCIObj );
 FLA_Obj_create_without_buffer( datatype, ldt, k, &tCOObj );
 FLA_Obj_attach_buffer( vCIBuff, 1, ldv, &vCIObj );
 FLA_Obj_attach_buffer( tauCIBuff, 1, k, &tauCIObj );
 FLA_Obj_attach_buffer( tCOBuff, 1, ldt, &tCOObj );

  //Call C function
  FLA_QR_blk_external( vCIObj, tauCIObj );
  larft_C( FLA_FORWARD, FLA_COLUMNWISE, vCIObj, tauCIObj, tCOObj );
  double diff =  computeError<T>( ldt, k, tCOBuff, tCPPOBuff );

  if(diff != 0.0)
  {
    printf( "larft(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "larft(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete vCPPIBuff ;
  delete tauCPPIBuff ;
  delete tCPPOBuff;
  FLA_Obj_free( &vCIObj );
  FLA_Obj_free( &tauCIObj );
  FLA_Obj_free( &tCOObj );
}

void larft_testall_variants(){
  larft_test<float>();
  larft_test<double>();
  larft_test<lapack_complex_float>();
  larft_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  larft_testall_variants();
}