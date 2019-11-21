
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

/*! @file liblame_potri.cc
 *  libflame_potri.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

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

    spotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    dpotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_COMPLEX:
  {
    lapack_complex_float *buff_A = ( lapack_complex_float * ) FLA_COMPLEX_PTR( A );

    cpotri_( &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                &info );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    lapack_complex_double *buff_A = ( lapack_complex_double * ) FLA_DOUBLE_COMPLEX_PTR( A );

    zpotri_( &blas_uplo,
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
  potri_C( FLA_LOWER_TRIANGULAR, aCIOObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, aCIOBuff, aCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "potri(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "potri(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete aCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
}

void potri_testall_variants(){
  potri_test<float>();
  potri_test<double>();
  potri_test<lapack_complex_float>();
  potri_test<lapack_complex_double>();
}

int main(int argc, char *argv[])
{
  potri_testall_variants();
}