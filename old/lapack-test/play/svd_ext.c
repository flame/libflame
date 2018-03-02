/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

typedef float testtype;
#define TESTTYPE FLA_FLOAT
#define REALTYPE FLA_FLOAT

//typedef scomplex testtype;
//#define TESTTYPE FLA_COMPLEX
//#define REALTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE
//#define REALTYPE FLA_DOUBLE

//typedef dcomplex testtype;
//#define TESTTYPE FLA_DOUBLE_COMPLEX
//#define REALTYPE FLA_DOUBLE

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Datatype realtype = REALTYPE;

  FLA_Svd_type jobu, jobv;
  FLA_Trans    transu, transv, transu_app, transv_app;
  FLA_Obj      A, s, U, V, A_copy, A_recovered;
  dim_t        m, n, min_m_n;
  FLA_Error    init_result; 
  double       residual = 0.0;
  int          test = 0;

  if ( argc == 3 || argc == 4 ) {
    m    = atoi(argv[1]);
    n    = atoi(argv[2]);
    if ( argc == 4 )
      test = atoi(argv[3]);
    min_m_n = min(m,n);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m n test\n", argv[0]);
    fprintf(stderr, "       m    : matrix length\n");
    fprintf(stderr, "       n    : matrix width\n");
    fprintf(stderr, "       test : test case \n");
    fprintf(stderr, "         0 - (all, all)\n");
    fprintf(stderr, "         1 - (mincopy, mincopy)\n");
    fprintf(stderr, "         2 - (mincopy, mincopy) with both conj transpose\n");
    fprintf(stderr, "         3 - (minover, mincopy)\n");
    fprintf(stderr, "         4 - (mincopy, minover)\n");
    fprintf(stderr, "         5 - (none, none)\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME SVD
  FLA_Obj_create( datatype, m, n, 0, 0, &A );
  FLA_Obj_create( realtype, min_m_n, 1, 0, 0, &s );

  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 
  FLA_Set_to_identity( A );
  {
    float* bA = FLA_Obj_buffer_at_view( A );
    bA[0] = -3.3675e-02; bA[3] = -4.7889e-01;  bA[6] = 1.2115e-01;
    bA[1] = -4.6349e-01;  bA[4] = 4.3787e-01;  bA[7] = 8.4413e-02;
    bA[2] = 6.6654e-02; bA[5] = -2.8100e-01; bA[8] = -9.0184e-01;
  }

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_recovered );

  switch ( test ) {
  case 0: 
    {
      jobu = FLA_SVD_VECTORS_ALL; transu = FLA_NO_TRANSPOSE; transu_app = FLA_NO_TRANSPOSE;
      jobv = FLA_SVD_VECTORS_ALL; transv = FLA_NO_TRANSPOSE; transv_app = FLA_CONJ_TRANSPOSE;

      FLA_Obj_create( datatype, m, m, 0, 0, &U  ); FLA_Set( FLA_ZERO, U ); 
      FLA_Obj_create( datatype, n, n, 0, 0, &V  ); FLA_Set( FLA_ZERO, V ); 
      break;
    }
  case 1:
    {
      jobu = FLA_SVD_VECTORS_MIN_COPY; transu = FLA_NO_TRANSPOSE; transu_app = FLA_NO_TRANSPOSE;
      jobv = FLA_SVD_VECTORS_MIN_COPY; transv = FLA_NO_TRANSPOSE; transv_app = FLA_CONJ_TRANSPOSE;

      FLA_Obj_create( datatype, m, min_m_n, 0, 0, &U  ); FLA_Set( FLA_ZERO, U ); 
      FLA_Obj_create( datatype, n, min_m_n, 0, 0, &V  ); FLA_Set( FLA_ZERO, V ); 
      break;
    }
  case 2:
    {
      jobu = FLA_SVD_VECTORS_MIN_COPY; transu = FLA_CONJ_TRANSPOSE; transu_app = FLA_CONJ_TRANSPOSE;
      jobv = FLA_SVD_VECTORS_MIN_COPY; transv = FLA_CONJ_TRANSPOSE; transv_app = FLA_NO_TRANSPOSE;

      FLA_Obj_create( datatype, min_m_n, m, 0, 0, &U  ); FLA_Set( FLA_ZERO, U ); 
      FLA_Obj_create( datatype, min_m_n, n, 0, 0, &V  ); FLA_Set( FLA_ZERO, V ); 
      break;
    }
  case 3:
    {
      jobu = FLA_SVD_VECTORS_MIN_OVERWRITE; transu = FLA_NO_TRANSPOSE; transu_app = FLA_NO_TRANSPOSE;
      jobv = FLA_SVD_VECTORS_MIN_COPY;      transv = FLA_NO_TRANSPOSE; transv_app = FLA_CONJ_TRANSPOSE;

      FLA_Obj_create( datatype, n, min_m_n, 0, 0, &V  ); FLA_Set( FLA_ZERO, V ); 
      break;
    }
  case 4:
    {
      jobu = FLA_SVD_VECTORS_MIN_COPY;      transu = FLA_NO_TRANSPOSE;   transu_app = FLA_NO_TRANSPOSE;
      jobv = FLA_SVD_VECTORS_MIN_OVERWRITE; transv = FLA_CONJ_TRANSPOSE; transv_app = FLA_NO_TRANSPOSE;

      FLA_Obj_create( datatype, m, min_m_n, 0, 0, &U  ); FLA_Set( FLA_ZERO, U ); 
      break;
    }
  case 5:
    {
      jobu = FLA_SVD_VECTORS_NONE;  transu = FLA_NO_TRANSPOSE;   transu_app = FLA_NO_TRANSPOSE;      
      jobv = FLA_SVD_VECTORS_NONE;  transu = FLA_NO_TRANSPOSE;   transu_app = FLA_CONJ_TRANSPOSE;
      // U and V are not created.
      break;
    }
  }
  FLA_Obj_fshow( stdout, " - s - ", s, "% 6.4e", "------");

  // Extended interface
  FLA_Svd_ext( jobu, transu, 
               jobv, transv,
               A, s, U, V );

  // U or V is A if A is overwritten; thus, properly partition U and V for checking
  {
    FLA_Obj W;
    if ( jobu == FLA_SVD_VECTORS_MIN_OVERWRITE ) {
      U = A;
      FLA_Part_1x2( U, &U, &W, min_m_n, FLA_LEFT ); 
    }
    if ( jobv == FLA_SVD_VECTORS_MIN_OVERWRITE ) {
      V = A;
      FLA_Part_2x1( V, &V, &W, min_m_n, FLA_TOP ); 
    }
  }

  FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
  FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
  switch( test ) {
  case 0:
  case 1:
  case 2:
    FLA_Obj_fshow( stdout, " - U - ", U, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - V - ", V, "% 6.4e", "------");
    break;
  case 3:
    FLA_Obj_fshow( stdout, " - A (U overwritten) - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - V - ", V, "% 6.4e", "------");
    break;
  case 4:
    FLA_Obj_fshow( stdout, " - U - ", U, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - A (V overwritten) - ", A, "% 6.4e", "------");
    break;
  }
  FLA_Obj_fshow( stdout, " - s - ", s, "% 6.4e", "------");

  // Recover the matrix
  if ( test < 5 ) {
    fprintf( stdout, "lapack2flame: %lu x %lu:\n ", m, n);
    // s is real, always no conjugate when it is used in apply_diag_matrix.
    if ( m >= n ) {
      FLA_Obj U1, U2;
      if ( transu_app == FLA_NO_TRANSPOSE ) {
        FLA_Part_1x2( U, &U1, &U2, n, FLA_LEFT );
        FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, U1 ); 
      } else {
        FLA_Part_2x1( U, &U1, &U2, n, FLA_TOP );
        FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, s, U1 );
      }
      FLA_Gemm_external( transu_app, transv_app,
                         FLA_ONE, U1, V, FLA_ZERO, A_recovered );
    } else {
      FLA_Obj V1, V2;
      if ( transv_app == FLA_CONJ_TRANSPOSE ) {
        FLA_Part_1x2( V, &V1, &V2, m, FLA_LEFT );
        FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, V1 );
      } else { 
        FLA_Part_2x1( V, &V1, &V2, m, FLA_TOP );
        FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, s, V1 ); 
      }
      FLA_Gemm_external( transu_app, transv_app,
                         FLA_ONE, U, V1, FLA_ZERO, A_recovered );
    }
    residual = FLA_Max_elemwise_diff( A_copy, A_recovered );
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");

    fprintf( stdout, "residual = %12.10e\n\n", residual );
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &s );

  if ( test != 3 && test != 5 ) 
    FLA_Obj_free( &U );
  if ( test != 4 && test != 5 )
    FLA_Obj_free( &V );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );

  FLA_Finalize_safe( init_result );     
}
