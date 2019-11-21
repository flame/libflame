
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

/*! @file liblame_stedc.cc
 *  libflame_stedc.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

template< typename T >
void stedc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  T *dCPPIOBuff, *eCPPIOBuff, *zCPPIOBuff;
  T *dCIOBuff, *eCIOBuff, *zCIOBuff;
  char jobz = 'V';
  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::stedc( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatype, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevdd_external( FLA_EVD_WITH_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<T>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<T>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<T>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stedc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stedc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );
}

template< typename Ta, typename Tb >
void stedc_test()
{
  int n = 512;
  srand (time(NULL));

  FLA_Init( );
  FLA_Obj dCIOObj, eCIOObj, zCIOObj;
  Tb *dCPPIOBuff, *eCPPIOBuff;
  Ta *zCPPIOBuff;
  Tb *dCIOBuff, *eCIOBuff;
  Ta *zCIOBuff;
  char jobz = 'V';
  int datatype = getDatatype<Ta>();
  int datatypeTb = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(dCPPIOBuff, dCIOBuff, n);
  allocate_init_buffer(eCPPIOBuff, eCIOBuff, (n-1));
  allocate_init_buffer(zCPPIOBuff, zCIOBuff, n*n);

  //Call CPP function
  libflame::stedc( LAPACK_COL_MAJOR, &jobz, &n, dCPPIOBuff, eCPPIOBuff, zCPPIOBuff, &n  );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatypeTb, n, 1, &dCIOObj );
  FLA_Obj_create_without_buffer( datatypeTb, (n-1), 1, &eCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &zCIOObj );
  FLA_Obj_attach_buffer( dCIOBuff, 1, n, &dCIOObj );
  FLA_Obj_attach_buffer( eCIOBuff, 1, (n-1), &eCIOObj );
  FLA_Obj_attach_buffer( zCIOBuff, 1, n, &zCIOObj );

  //Call C function
  FLA_Tevdd_external( FLA_EVD_WITH_VECTORS, dCIOObj, eCIOObj, zCIOObj);
  double diff =  computeError<Tb>( n, 1, dCIOBuff, dCPPIOBuff );
  diff +=  computeError<Tb>( (n-1), 1, eCIOBuff, eCPPIOBuff );
  diff +=  computeError<Ta>( n, n, zCIOBuff, zCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "stedc(): %s TEST FAIL\n" , __PRETTY_FUNCTION__);
  }else{
    printf( "stedc(): %s TEST PASS\n" , __PRETTY_FUNCTION__);
  }

  //Free up the buffers
  delete dCPPIOBuff ;
  delete eCPPIOBuff ;
  delete zCPPIOBuff ;
  FLA_Obj_free( &dCIOObj );
  FLA_Obj_free( &eCIOObj );
  FLA_Obj_free( &zCIOObj );
}

void stedc_testall_variants(){
  stedc_test<float>();
  stedc_test<double>();
  stedc_test<lapack_complex_float, float>();
  stedc_test<lapack_complex_double, double >();
}

int main(int argc, char *argv[])
{
  stedc_testall_variants();
}