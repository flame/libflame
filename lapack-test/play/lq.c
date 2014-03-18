#include "FLAME.h"
#define EPS 1.e-10

//typedef float testtype;
//#define TESTTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE

typedef dcomplex testtype;
#define TESTTYPE FLA_DOUBLE_COMPLEX

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Obj      A, Ak, T, Tk, D, Dk, A_copy, A_recovered, L, Q, Qk, W, x, y, z;
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

  // FLAME LQ^H setup
  FLA_Obj_create( datatype, m, n, 0, 0, &A );
  FLA_LQ_UT_create_T( A, &T );
  
  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );
  FLA_Copy( A, A_copy );

  // LQ test ( A = L Q^H )
  FLA_LQ_UT( A, T );

  // Create Q (identity), L (A_copy)
  FLA_Obj_create( datatype, m, n, 0, 0, &Q  ); FLA_Set_to_identity( Q  );
  FLA_Obj_create( datatype, m, m, 0, 0, &D  );

  FLA_Obj_create( datatype, k, n, 0, 0, &Qk ); FLA_Set_to_identity( Qk );
  FLA_Obj_create( datatype, k, k, 0, 0, &Dk  );

  
  FLA_Obj_create( datatype, m, m, 0, 0, &L );

  // Q^H := I H_{0}^H ... H_{k-1}^H
  if ( use_form_q ) {
    FLA_LQ_UT_form_Q( A, T, Q );
  } else {
    FLA_Apply_Q_UT_create_workspace_side( FLA_RIGHT, T, Q, &W );
    FLA_Apply_Q_UT( FLA_RIGHT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                    A, T, W, Q );
    FLA_Obj_free( &W );
  }

  // D := Q^T Q
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, 
                     FLA_ONE, Q, Q, FLA_ZERO, D ); 

  // Qk := I H0 ... Hk
  FLA_Part_1x2( T, &Tk, &W, k, FLA_LEFT );
  FLA_Part_2x1( A, &Ak, &W, k, FLA_TOP );

  if ( use_form_q ) {
    // Overwrite the result to test FLAME API
    FLA_Set( FLA_ZERO, Qk );
    FLA_Copy( Ak, Qk );
    FLA_LQ_UT_form_Q( Ak, Tk, Qk );
  } else {
    FLA_Apply_Q_UT_create_workspace( Tk, Qk, &W );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                    Ak, Tk, W, Qk );
    FLA_Obj_free( &W );
  }

  // Dk := Qk^T Qk
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, 
                     FLA_ONE, Qk, Qk, FLA_ZERO, Dk ); 

  // L := A (Q^H)^H
  if ( use_form_q ) {
    // Note that the formed Q is actually Q^H; transb should be carefully assigned.
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, 
                       FLA_ONE, A_copy, Q, FLA_ZERO, L );
  } else {
    FLA_Apply_Q_UT_create_workspace( T, L, &W );
    FLA_Apply_Q_UT( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                    A, T, W, L );
    FLA_Obj_free( &W );
  }
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, L, Q, FLA_ZERO, A_recovered ); 

  // Create vectors for testing
  FLA_Obj_create( datatype, n, 1, 0, 0, &x ); FLA_Set( FLA_ZERO, x );
  FLA_Obj_create( datatype, m, 1, 0, 0, &y ); FLA_Set( FLA_ZERO, y );
  FLA_Obj_create( datatype, m, 1, 0, 0, &z ); FLA_Set( FLA_ZERO, z );

  // x is given
  FLA_Set( FLA_ONE, x );

  // y := Ax
  FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_copy, x, FLA_ZERO, y );
  
  // z := L (Q^H) x , libflame
  FLA_Apply_Q_UT_create_workspace( T, x, &W );
  FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                  A, T, W, x );
  FLA_Obj_free( &W );

  if ( m < n )
    FLA_Part_2x1( x, &x, &W, m, FLA_TOP );
  else 
    FLA_Part_1x2( L, &L, &W, n, FLA_LEFT );    

  FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, L, x, FLA_ZERO, z );
  
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
    FLA_Obj_fshow( stdout, " - L - ", L, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu, %lu: ", m, n, k);
    fprintf( stdout, "| A - A_recovered | = %12.10e, | Ax - y | = %12.10e\n\n",
             residual_A, residual_Axy ) ;
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &T );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );
  FLA_Obj_free( &L );
  FLA_Obj_free( &Q );
  FLA_Obj_free( &Qk );
  FLA_Obj_free( &D );
  FLA_Obj_free( &Dk );

  FLA_Obj_free( &x );
  FLA_Obj_free( &y );
  FLA_Obj_free( &z );
                   
  FLA_Finalize_safe( init_result );     
}
