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
  FLA_Obj      
    A, TU, TV, 
    A_copy, A_recovered,
    G, H, U, V, Vb, d, e, 
    DU, DV;

  FLA_Uplo     uplo;
  dim_t        m, n, min_m_n;
  FLA_Error    init_result; 

  double       residual_A = 0.0;

  if ( argc == 3 ) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    min_m_n = min(m,n);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m n k\n", argv[0]);
    fprintf(stderr, "       m : matrix length\n");
    fprintf(stderr, "       n : matrix width\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 || n == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME Bidiag setup
  FLA_Obj_create( datatype, m, n, 0, 0, &A );

  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 

  uplo = ( m >= n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );

  // Bidiag test
  FLA_Bidiag_UT_create_T( A, &TU, &TV );
  FLA_Bidiag_UT( A, TU, TV );

  // Orthonomal basis U, V. 
  FLA_Obj_create( datatype, m, min_m_n, 0, 0, &U ); 
  FLA_Obj_create( datatype, n, min_m_n, 0, 0, &V ); 

  FLA_Set( FLA_ZERO, U );
  FLA_Set( FLA_ZERO, V );

  FLA_Bidiag_UT_form_U( A, TU, U );
  FLA_Bidiag_UT_form_V( A, TV, V ); 

  FLA_Obj_free( &TU );
  FLA_Obj_free( &TV );
  
  if ( FLA_Obj_is_complex( A ) ){
    FLA_Obj rL, rR;
    FLA_Obj U1, U2;
    FLA_Obj V1, V2;
    
    FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &rL );
    FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &rR );

    FLA_Bidiag_UT_realify( A, rL, rR );

    FLA_Obj_fshow( stdout, " - rL - ", rL, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - rR - ", rR, "% 6.4e", "------");

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

  // Create workspace to store given's scalars.
  FLA_Bsvd_create_workspace( d, &G, &H );

  // BSVD
  FLA_Bsvd( uplo, d, e, G, H,
            FLA_SVD_VECTORS_ALL, U,
            FLA_SVD_VECTORS_ALL, V );

  FLA_Obj_free( &G );
  FLA_Obj_free( &H );

  if ( min_m_n > 1 )
    FLA_Obj_free( &e );

  // U^H U
  FLA_Obj_create( datatype, min_m_n, min_m_n, 0,0, &DU );
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, U, U, FLA_ZERO, DU );

  // V^H V
  FLA_Obj_create( datatype, min_m_n, min_m_n, 0,0, &DV );


 
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, V, V, FLA_ZERO, DV );
  
  // Recover the matrix

  // Vb := diag(d) (V^H)
  FLA_Obj_create_copy_of( FLA_CONJ_TRANSPOSE, V, &Vb );
  FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, d, Vb );

  // A := U Vb
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, U, Vb, FLA_ZERO, A_recovered );

  residual_A    = FLA_Max_elemwise_diff( A_copy, A_recovered );

  if (1) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - d - ", d, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - U - ", U, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - V - ", V, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Vb - ", Vb, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - U'U - ", DU,  "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - VV' - ", DV,  "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu: ", m, n);
    fprintf( stdout, "recovery A = %12.10e\n\n", residual_A ) ;
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &d );

  FLA_Obj_free( &U );
  FLA_Obj_free( &V );
  FLA_Obj_free( &Vb );

  FLA_Obj_free( &DU );
  FLA_Obj_free( &DV );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );


  FLA_Finalize_safe( init_result );     
}
