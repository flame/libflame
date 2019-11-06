
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

/*! @file liblame_trsyl.cc
 *  libflame_trsyl.cc Test application to validate CPP template interface
 *  */
 
#include "libflame_test.hh"

extern TLS_CLASS_SPEC FLA_Obj FLA_ONE;
template< typename T >
void trsyl_test()
{
  int m = 64;
  int n = 64;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  FLA_Trans transARef = FLA_NO_TRANSPOSE ; FLA_Trans transBRef = FLA_NO_TRANSPOSE;
  FLA_Obj aCIOObj,  bCIOObj, cCIOObj, scaleRefObj ;
  T *aCPPIOBuff, *bCPPIOBuff, *cCPPIOBuff ;
  T *aCIOBuff, *bCIOBuff, *cCIOBuff, scaleRef ;
  T scale;
  int isgn = 1;

  int datatype = getDatatype<T>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);

  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bCPPIOBuff, &n, cCPPIOBuff, &m, &scale );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, n, &bCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( &scaleRef, 1, 1, &scaleRefObj );

  //Call C function
  FLA_Sylv_unb_external( transARef, transBRef, FLA_ONE, aCIOObj,  bCIOObj, cCIOObj, scaleRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<T>( m, m, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): Failure Diff = %E\n", diff);
  }else{
   printf( "trsyl(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bCPPIOBuff;
  delete cCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
  FLA_Obj_free( &cCIOObj );
}


template< typename Ta, typename Tb  >
void trsyl_test()
{
  int m = 64;
  int n = 64;
  srand (time(NULL));

  FLA_Init( );
  char transa = 'N' ;char transb = 'N';
  FLA_Trans transARef = FLA_NO_TRANSPOSE ; FLA_Trans transBRef = FLA_NO_TRANSPOSE;
  FLA_Obj aCIOObj,  bCIOObj, cCIOObj, scaleRefObj ;
  Tb *aCPPIOBuff, *bCPPIOBuff, *cCPPIOBuff ;
  Tb *aCIOBuff, *bCIOBuff, *cCIOBuff, scaleRef ;
  Tb scale;
  int isgn = 1;

  int datatype = getDatatype<Tb>();

  //Allocate and initialize buffers for C and CPP functions with random values
  allocate_init_buffer(aCPPIOBuff, aCIOBuff, m*m);
  allocate_init_buffer(bCPPIOBuff, bCIOBuff, n*n);
  allocate_init_buffer(cCPPIOBuff, cCIOBuff, m*n);

  //Call CPP function
  libflame::trsyl( LAPACK_COL_MAJOR, &transa, &transb, &isgn, &m, &n, aCPPIOBuff, &m, bCPPIOBuff, &n, cCPPIOBuff, &m, &scale );

  //Allocate Object for C function and copy already allocated and filled buffer
  FLA_Obj_create_without_buffer( datatype, m, m, &aCIOObj );
  FLA_Obj_attach_buffer( aCIOBuff, 1, m, &aCIOObj );
  FLA_Obj_create_without_buffer( datatype, n, n, &bCIOObj );
  FLA_Obj_attach_buffer( bCIOBuff, 1, n, &bCIOObj );
  FLA_Obj_create_without_buffer( datatype, m, n, &cCIOObj );
  FLA_Obj_attach_buffer( cCIOBuff, 1, m, &cCIOObj );
  FLA_Obj_create_without_buffer( datatype, 1, 1, &scaleRefObj );
  FLA_Obj_attach_buffer( &scaleRef, 1, 1, &scaleRefObj );

  //Call C function
  FLA_Sylv_unb_external( transARef, transBRef, FLA_ONE, aCIOObj,  bCIOObj, cCIOObj, scaleRefObj );

  //Compute Difference in C and CPP buffer
  double diff =  computeError<Tb>( m, m, cCIOBuff, cCPPIOBuff );

  if(diff != 0.0)
  {
    printf( "trsyl(): Failure Diff = %E\n", diff);
  }else{
   printf( "trsyl(): Success\n");
  }

  //Free up the buffers
  delete aCPPIOBuff;
  delete bCPPIOBuff;
  delete cCPPIOBuff;
  FLA_Obj_free( &aCIOObj );
  FLA_Obj_free( &bCIOObj );
  FLA_Obj_free( &cCIOObj );
}

void trsyl_testall_variants(){
  trsyl_test<float>();
  trsyl_test<double>();
  trsyl_test<lapack_complex_float, float>();
  trsyl_test<lapack_complex_double, double>();
}

int main(int argc, char *argv[])
{
  trsyl_testall_variants();
}