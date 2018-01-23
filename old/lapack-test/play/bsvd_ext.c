/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"
#define EPS 1.e-10

//typedef float testtype;
//#define TESTTYPE FLA_FLOAT
//#define REALTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE
//#define REALTYPE FLA_DOUBLE

typedef scomplex testtype;
#define TESTTYPE FLA_COMPLEX
#define REALTYPE FLA_FLOAT

//typedef dcomplex testtype;
//#define TESTTYPE FLA_DOUBLE_COMPLEX
//#define REALTYPE FLA_DOUBLE

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Datatype realtype = REALTYPE;
  FLA_Obj      A, d, e, U, V, C, D, W, A_copy, C_copy;
  FLA_Uplo     uplo;
  dim_t        m, n, min_m_n;
  FLA_Error    init_result; 
  int          test = 2;
  double       residual = 0.0;

  if ( argc == 4 ) {
    m    = atoi(argv[1]);
    n    = atoi(argv[2]);
    test = atoi(argv[3]);
    min_m_n = min(m,n);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m n test\n", argv[0]);
    fprintf(stderr, "       m    : matrix length\n");
    fprintf(stderr, "       n    : matrix width\n");
    fprintf(stderr, "       test : test case\n");
    fprintf(stderr, "         0 - (all, all)\n");
    fprintf(stderr, "         1 - (none, all)\n");
    fprintf(stderr, "         2 - (none, none)\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME Bidiag setup
  FLA_Obj_create( datatype, m, n, 0, 0, &A );       // LHS
  FLA_Obj_create( datatype, m, 1, 0, 0, &D );       // RHS
  FLA_Obj_create( datatype, n, 1, 0, 0, &C );       // Solution vector
  FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &W ); // Work space

  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );

  FLA_Random_matrix( C );
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, C, &C_copy );

  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, A, C, FLA_ZERO, D );

  uplo = ( m >= n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );

  // Orthonomal basis U, V. 
  FLA_Obj_create( datatype, m, min_m_n, 0, 0, &U ); 
  FLA_Obj_create( datatype, n, min_m_n, 0, 0, &V ); 

  // Bidiag test
  {
    FLA_Obj TU, TV;
    FLA_Bidiag_UT_create_T( A, &TU, &TV );
    FLA_Bidiag_UT( A, TU, TV );

    FLA_Set( FLA_ZERO, U );
    FLA_Set( FLA_ZERO, V );

    FLA_Bidiag_UT_form_U( A, TU, U );
    FLA_Bidiag_UT_form_V( A, TV, V ); 

    FLA_Obj_free( &TU );
    FLA_Obj_free( &TV );
  }

  // Realify U and V
  if ( FLA_Obj_is_complex( A ) ){
    FLA_Obj rL, rR;
    FLA_Obj U1, U2;
    FLA_Obj V1, V2;
    
    FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &rL );
    FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &rR );

    FLA_Bidiag_UT_realify( A, rL, rR );

    FLA_Part_1x2( U,   &U1, &U2,   min_m_n, FLA_LEFT );
    FLA_Part_1x2( V,   &V1, &V2,   min_m_n, FLA_LEFT );
    
    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE,    rL, U1 );
    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, V1 );

    FLA_Obj_free( &rL );
    FLA_Obj_free( &rR );
  }

  // Extract diagonals
  FLA_Obj_create( realtype, min_m_n, 1, 0, 0, &d );  
  if ( min_m_n > 1 ) 
    FLA_Obj_create( realtype, min_m_n-1 , 1, 0, 0, &e );  

  FLA_Bidiag_UT_extract_real_diagonals( A, d, e );

  FLA_Obj_fshow( stdout, " - U Bidiag- ", U, "% 6.4e", "------");

  // BSVD
  {
    FLA_Obj G, H;
    FLA_Svd_type jobu, jobv;
    FLA_Bool apply_Uh2C;

    switch( test ) {
    case 0: {
      jobu = FLA_SVD_VECTORS_ALL;
      jobv = FLA_SVD_VECTORS_ALL;
      apply_Uh2C = FALSE;
      break;
    }
    case 1: {
      jobu = FLA_SVD_VECTORS_NONE;
      jobv = FLA_SVD_VECTORS_ALL;
      apply_Uh2C = TRUE;
      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, U, D, FLA_ZERO, W );
      break;
    }
    case 2: {
      jobu = FLA_SVD_VECTORS_NONE;
      jobv = FLA_SVD_VECTORS_NONE;
      apply_Uh2C = FALSE;
      break;
    }
    }
    
    // Create workspace to store given's scalars.
    FLA_Bsvd_create_workspace( d, &G, &H );
    FLA_Bsvd_ext( uplo, d, e, G, H,
                  jobu, U,
                  jobv, V, 
                  apply_Uh2C, W );
    FLA_Obj_free( &G );
    FLA_Obj_free( &H );
  }

  if ( min_m_n > 1 )
    FLA_Obj_free( &e );

  if (1) {
    //FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - d - ", d, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - U - ", U, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - V - ", V, "% 6.4e", "------");
  }

  // Equating Ax = b
  {
    switch( test ) {
    case 0: {
      FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, U, D, FLA_ZERO, W );
      FLA_Obj_fshow( stdout, " - W - ", W, "% 6.4e", "------");
      FLA_Invert( FLA_NO_CONJUGATE, d );
      FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, d, W );
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, V, W, FLA_ZERO, C );
      break;
    }
    case 1: {
      FLA_Obj_fshow( stdout, " - W - ", W, "% 6.4e", "------");
      FLA_Invert( FLA_NO_CONJUGATE, d );
      FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, d, W );
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, V, W, FLA_ZERO, C );
      break;
    }
    case 2: {
      // Do nothing
      break;
    }
    }
  }
  residual = FLA_Max_elemwise_diff( C, C_copy );
  FLA_Obj_free( &W );

  if (1) {
    FLA_Obj_fshow( stdout, " - Solution - ", C, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu: ", m, n);
    fprintf( stdout, "residual = %12.10e\n\n", residual ) ;
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &D );
  FLA_Obj_free( &C );
  FLA_Obj_free( &d );

  FLA_Obj_free( &U );
  FLA_Obj_free( &V );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &C_copy );


  FLA_Finalize_safe( init_result );     
}
