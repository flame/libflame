/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"
#define EPS 1.e-10

#define TESTTYPE FLA_DOUBLE
#define REALTYPE FLA_DOUBLE

int main( int argc, char** argv ) {
  FLA_Obj      A, T, p, w, D, A_copy, A_recovered, R, Q, W;
  dim_t        m, n;
  FLA_Error    init_result; 
  double       residual_A, residual_R;

  if ( argc == 3 ) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m n\n", argv[0]);
    fprintf(stderr, "       m : matrix length\n");
    fprintf(stderr, "       n : matrix width\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME QR piv setup
  FLA_Obj_create( TESTTYPE, m,n, 0,0, &A );
  FLA_Obj_create( REALTYPE, n,1, 0,0, &w );
  FLA_Obj_create( FLA_INT,  n,1, 0,0, &p );
  FLA_QR_UT_create_T( A, &T );
  
  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );
  FLA_Copy( A, A_copy );

  // QR test; AP^T = QR
  FLA_QR_UT_piv( A, T, w, p );

  // Q := H0 ... H_all 
  FLA_Obj_create( TESTTYPE, m, m, 0,0, &Q  ); 
  FLA_Set_to_identity( Q );

  FLA_Apply_Q_UT_create_workspace( T, Q, &W );
  FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                  A, T, W, Q );
  FLA_Obj_free( &W );

  // D := Q^T Q
  FLA_Obj_create( TESTTYPE, m, m, 0,0, &D  );
  FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, Q, Q, FLA_ZERO, D ); 

  // R := Q^T AP^T
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A_copy, &R );
  FLA_Copy( A_copy, R );

  FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE, p, R );
  FLA_Apply_Q_UT_create_workspace( T, R, &W );
  FLA_Apply_Q_UT( FLA_LEFT, FLA_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                  A, T, W, R );
  FLA_Obj_free( &W );

  // Recover A
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, Q, R, FLA_ZERO, A_recovered ); 

  // Pivot
  FLA_Apply_pivots( FLA_RIGHT, FLA_NO_TRANSPOSE, p, A_recovered );
  
  // Comapre (A_copy, A_recovered)
  residual_A      = FLA_Max_elemwise_diff( A_copy, A_recovered );

  if (1 || residual_A > EPS ) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - T - ", T, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - p - ", p, "%d", "------");
    FLA_Obj_fshow( stdout, " - w - ", w, "%e", "------");
    FLA_Obj_fshow( stdout, " - Q - ", Q, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - D = Q^T Q - ", D, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - R - ", R, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
  }
  FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
  residual_R      = FLA_Max_elemwise_diff( A, R );

  fprintf( stdout, "lapack2flame: %lu x %lu: \n", m, n);
  fprintf( stdout, "recovery A = %12.10e\n\n", residual_A);
  fprintf( stdout, "recovery R = %12.10e\n\n", residual_R);
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &T );
  FLA_Obj_free( &p );
  FLA_Obj_free( &w );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );
  FLA_Obj_free( &R );
  FLA_Obj_free( &Q );
  FLA_Obj_free( &D );
                   
  FLA_Finalize_safe( init_result );     
}
