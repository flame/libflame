#include "FLAME.h"
#define EPS 1.e-10

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
  FLA_Obj      A, s, U, V, DU, DV, A_copy, A_recovered, norm;
  dim_t        m, n, min_m_n;
  FLA_Error    init_result; 
  double       residual_A = 0.0;
  int          test = 0;

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
    fprintf(stderr, "       test : test case 0 - all, all, 1 - mincopy, mincopy\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME SVD
  FLA_Obj_create( datatype, m, n, 0, 0, &A );
  FLA_Obj_create( realtype, min_m_n, 1, 0, 0, &s );
  FLA_Obj_create( realtype, 1, 1, 0, 0, &norm );

  switch ( test ) {
  case 0: 
    {
      jobu = FLA_SVD_VECTORS_ALL;
      jobv = FLA_SVD_VECTORS_ALL;

      FLA_Obj_create( datatype, m, m, 0, 0, &U  );
      FLA_Obj_create( datatype, m, m, 0, 0, &DU );

      FLA_Obj_create( datatype, n, n, 0, 0, &V  );
      FLA_Obj_create( datatype, n, n, 0, 0, &DV );
      break;
    }
  case 1:
    {
      jobu = FLA_SVD_VECTORS_MIN_COPY;
      jobv = FLA_SVD_VECTORS_MIN_COPY;

      FLA_Obj_create( datatype, m, min_m_n, 0, 0, &U  );
      FLA_Obj_create( datatype, min_m_n, min_m_n, 0, 0, &DU );

      FLA_Obj_create( datatype, n, min_m_n, 0, 0, &V  );
      FLA_Obj_create( datatype, min_m_n, min_m_n, 0, 0, &DV );
      break;
    }
  }

  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );

  FLA_Set( FLA_ZERO, U );  FLA_Set_to_identity( DU ); 
  FLA_Set( FLA_ZERO, V );  FLA_Set_to_identity( DV ); 

  // SVD test
  FLA_Svd( jobu, jobv, A, s, U, V );

  if (1) {
    FLA_Obj_fshow( stdout, " - U - ", U, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - V - ", V, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - s - ", s, "% 6.4e", "------");
  }

  // Test orthonomal basis U^H U = I
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_MINUS_ONE, U, U, FLA_ONE, DU );
  FLA_Norm_frob( DU, norm );
  FLA_Obj_fshow(stdout, " - norm U - ", norm, "%6.4e", "--");

  // Test orthonomal basis V^H V = I
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_MINUS_ONE, V, V, FLA_ONE, DV );
  FLA_Norm_frob( DV, norm );
  FLA_Obj_fshow(stdout, " - norm V - ", norm, "%6.4e", "--");
  
  // Recover the matrix
  if ( m >= n ) {
    FLA_Obj UL, UR;
    FLA_Part_1x2( U, &UL, &UR, n, FLA_LEFT );
    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, UL, V, FLA_ZERO, A_recovered );
  } else {
    FLA_Obj VL, VR;
    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, U );
    FLA_Part_1x2( V, &VL, &VR, m, FLA_LEFT );
    FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, U, VL, FLA_ZERO, A_recovered );
  }
  residual_A    = FLA_Max_elemwise_diff( A_copy, A_recovered );

  if (1) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    //FLA_Obj_fshow( stdout, " - s - ", s, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu: ", m, n);
    fprintf( stdout, "recovery A = %12.10e\n\n", residual_A );
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &s );
  FLA_Obj_free( &U );
  FLA_Obj_free( &V );

  FLA_Obj_free( &DU );
  FLA_Obj_free( &DV );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );
  FLA_Obj_free( &norm );

  FLA_Finalize_safe( init_result );     
}
