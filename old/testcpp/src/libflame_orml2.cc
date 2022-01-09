
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

/*! @file liblame_orml2.cc
 *  libflame_orml2.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

FLA_Error unml2_C( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          cs_A;
  int          cs_B;
  int          k_t;
  int          lwork;
  char         blas_side;
  char         blas_trans;
  FLA_Obj      work_obj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Apply_Q_check( side, trans, storev, A, t, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  cs_B     = FLA_Obj_col_stride( B );

  k_t      = FLA_Obj_vector_dim( t );

  FLA_Param_map_flame_to_netlib_side( side, &blas_side );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );


  if ( side == FLA_LEFT )
    lwork  = n_B * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  else
    lwork  = m_B * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );

  switch( datatype ){
    case FLA_FLOAT:
    {
      float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_t    = ( float * ) FLA_FLOAT_PTR( t );
      float *buff_B    = ( float * ) FLA_FLOAT_PTR( B );
      float *buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );

       sormlq_( &blas_side,
                   &blas_trans,
                   &m_B,
                   &n_B,
                   &k_t,
                   buff_A, &cs_A,
                   buff_t,
                   buff_B, &cs_B,
                   buff_work, &lwork,
                   &info );

      break;
    }
    case FLA_DOUBLE:
    {
      double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
      double *buff_B    = ( double * ) FLA_DOUBLE_PTR( B );
      double *buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

     dormlq_( &blas_side,
                 &blas_trans,
                 &m_B,
                 &n_B,
                 &k_t,
                 buff_A, &cs_A,
                 buff_t,
                 buff_B, &cs_B,
                 buff_work, &lwork,
                 &info );

      break;
    }

    case FLA_COMPLEX:
    {
      lapack_complex_float *buff_A    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );
      lapack_complex_float *buff_t    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( t );
      lapack_complex_float *buff_B    = ( lapack_complex_float * ) FLA_COMPLEX_PTR( B );
      lapack_complex_float *buff_work = ( lapack_complex_float * ) FLA_COMPLEX_PTR( work_obj );

      cunmlq_( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
                  &info );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      lapack_complex_double *buff_A    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );
      lapack_complex_double *buff_t    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( t );
      lapack_complex_double *buff_B    = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( B );
      lapack_complex_double *buff_work = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

      zunmlq_( &blas_side,
                  &blas_trans,
                  &m_B,
                  &n_B,
                  &k_t,
                  buff_A, &cs_A,
                  buff_t,
                  buff_B, &cs_B,
                  buff_work, &lwork,
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

template< typename T >
void orml2_test()
{
  int m = 512;
  int n = 256;
  int k = 128;
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
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, k*n);
  tauCPPOBuff =  new T [k];
  tauCOBuff =  new T [k];
  cCIOBuff =  new T [m*n];
  cCPPIOBuff =  new T [m*n];
  for (int i = 0; i < m*n; i++)
  {
    cCPPIOBuff[i] = 0;
    cCIOBuff[i] = 0;
  }

  //Call CPP function
  libflame::gelqf( LAPACK_COL_MAJOR, &k, &n, aCPPIOBuff, &k, tauCPPOBuff );
  libflame::orml2( LAPACK_COL_MAJOR, &side, &trans, &m, &n, &k, aCPPIOBuff, &k, tauCPPOBuff, cCPPIOBuff, &m );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, k, n, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, k, 1, &tauCOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, k, &aCIOObj );
  FLA_Obj_attach_buffer( tauCOBuff, 1, k, &tauCOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );

  //Call C function
  FLA_LQ_blk_external( aCIOObj, tauCOObj );
  unml2_C( FLA_RIGHT, FLA_TRANSPOSE, FLA_ROWWISE, aCIOObj, tauCOObj, cCIOObj );

  double diff =  computeError<T>( k, n, aCIOBuff, aCPPIOBuff );
  diff +=  computeError<T>( 1, k, tauCOBuff, tauCPPOBuff );
  diff +=  computeError<T>( m, n, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "orml2(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "orml2(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete[] aCPPIOBuff ;
  delete[] tauCPPOBuff ;
  delete[] aCIOBuff ;
  delete[] tauCOBuff ;
  FLA_Obj_free_without_buffer( &aCIOObj );
  FLA_Obj_free_without_buffer( &tauCOObj );

}

void orml2_testall_variants(){
  orml2_test<float>();
  orml2_test<double>();
}

int main(int argc, char *argv[])
{
  orml2_testall_variants();
}