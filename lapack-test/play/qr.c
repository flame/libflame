#include "FLAME.h"
#define EPS 1.e-10

typedef double testtype;
#define TESTTYPE FLA_DOUBLE

//typedef float testtype;
//#define TESTTYPE FLA_FLOAT



int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Obj      A, Ak, T, Tk, D, Dk, A_copy, A_recovered, R, Q, Qk, W, x, y, z;
  dim_t        m, n, k;
  dim_t        min_m_n;
  FLA_Error    init_result; 
  double       residual_A, residual_Axy;
  int          use_form_q = 1;
  if ( argc == 4 ) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    k = atoi(argv[3]);
    min_m_n = min(m,n);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m n k\n", argv[0]);
    fprintf(stderr, "       m : matrix length\n");
    fprintf(stderr, "       n : matrix width\n");
    fprintf(stderr, "       k : number of house holder vectors applied for testing\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME QR setup
  FLA_Obj_create( datatype, m,n, 0,0, &A );
  FLA_QR_UT_create_T( A, &T );
  
  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 
  //FLA_Set_to_identity( A ); 

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );
  FLA_Copy( A, A_copy );

  // QR test
  FLA_QR_UT( A, T );

  // Create Q (identity), R (A_copy)
  FLA_Obj_create( datatype, m, m, 0,0, &Q  ); FLA_Set_to_identity( Q  );
  FLA_Obj_create( datatype, m, m, 0,0, &D  );

  FLA_Obj_create( datatype, m, k, 0,0, &Qk ); FLA_Set_to_identity( Qk );
  FLA_Obj_create( datatype, k, k, 0,0, &Dk  );

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_copy, &R );

  // Q := H0 ... H_all 
  if ( use_form_q ) {
    FLA_QR_UT_form_Q( A, T, Q );   
  } else {
    FLA_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, Q, &W );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                    A, T, W, Q );
    FLA_Obj_free( &W );
  }

  // D := Q^T Q
  FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, Q, Q, FLA_ZERO, D ); 

  // Qk := H0 ... Hk
  FLA_Part_1x2( T, &Tk, &W, k, FLA_LEFT );
  FLA_Part_1x2( A, &Ak, &W, k, FLA_LEFT );
  FLA_Apply_Q_UT_create_workspace_side( FLA_LEFT,  Tk, Qk, &W );
  FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                  Ak, Tk, W, Qk );
  FLA_Obj_free( &W );

  // Overwrite the result to test FLAME API
  if ( use_form_q ) {
    FLA_Set( FLA_ZERO, Qk );
    FLA_Copy( Ak, Qk );
    FLA_QR_UT_form_Q( Qk, Tk, Qk );   
  }

  // Dk := Qk^T Qk
  FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, Qk, Qk, FLA_ZERO, Dk ); 

  // R := Q^T A
  FLA_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, R, &W );
  if ( use_form_q ) {
    FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, Q, A_copy, FLA_ZERO, R );
  } else {
    FLA_Apply_Q_UT( FLA_LEFT, FLA_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                    A, T, W, R );
  }
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, Q, R, FLA_ZERO, A_recovered ); 

  // Create vectors for testing
  FLA_Obj_create( FLA_Obj_datatype( A ), n, 1, 0, 0, &x ); FLA_Set( FLA_ZERO, x );
  FLA_Obj_create( FLA_Obj_datatype( A ), m, 1, 0, 0, &y ); FLA_Set( FLA_ZERO, y );
  FLA_Obj_create( FLA_Obj_datatype( A ), m, 1, 0, 0, &z ); FLA_Set( FLA_ZERO, z );

  // x is given
  FLA_Set( FLA_ONE, x );

  // y := Ax
  FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_copy, x, FLA_ZERO, y );
  
  // z := QRx , libflame
  FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, R, x, FLA_ZERO, z );
  FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                  A, T, W, z );
  
  // Comapre (A_copy, A_recovered), (y,z) and (y,w)
  residual_A    = FLA_Max_elemwise_diff( A_copy, A_recovered );
  residual_Axy  = FLA_Max_elemwise_diff( y, z );

  if ( 1 || residual_A > EPS || residual_Axy > EPS ) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - T - ", T, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Q - ", Q, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - D = Q^T Q - ", D, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Qk - ", Qk, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Dk = Qk^T Qk - ", Dk, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - R - ", R, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu, %lu: ", m, n, k);
    fprintf( stdout, "| A - A_recovered | = %12.10e, | Ax - y | = %12.10e\n\n",
             residual_A, residual_Axy ) ;
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &T );


  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );
  FLA_Obj_free( &R );
  FLA_Obj_free( &Q );
  FLA_Obj_free( &Qk );
  FLA_Obj_free( &D );
  FLA_Obj_free( &Dk );
  FLA_Obj_free( &W );
  FLA_Obj_free( &x );
  FLA_Obj_free( &y );
  FLA_Obj_free( &z );
                   
  FLA_Finalize_safe( init_result );     
}
