#include "FLAME.h"
#define EPS 1.e-10

typedef float testtype;
#define TESTTYPE FLA_FLOAT
#define REALTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE
//#define REALTYPE FLA_DOUBLE

//typedef scomplex testtype;
//#define TESTTYPE FLA_COMPLEX
//#define REALTYPE FLA_FLOAT

//typedef dcomplex testtype;
//#define TESTTYPE FLA_DOUBLE_COMPLEX
//#define REALTYPE FLA_DOUBLE

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Datatype realtype = REALTYPE;
  FLA_Obj      
    A, TU, TV, 
    A_copy, A_recovered,
    U, V, Vb, B, Be, d, e, 
    DU, DV;

  FLA_Obj     
    ATL, ATR,
    ABL, ABR, Ae;

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
    fprintf(stderr, "Usage: %s m n\n", argv[0]);
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
  FLA_Bidiag_UT_create_T( A, &TU, &TV );

  // Rand A and create A_copy.
  FLA_Random_matrix( A ); 
  //FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
  //FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );


  // An identity matrix is a trivial case but is more tricky as 
  // its column and row are zeros (can easily incur 'divide by zero'
  // in house holder operations (e.g., re-arranged bidiagonalization 
  // process is not safe to do so.
  //FLA_Set_to_identity( A ); 

  if ( m >= n ) {
    uplo = FLA_UPPER_TRIANGULAR;
    FLA_Part_2x2( A, &ATL, &ATR,
                     &ABL, &ABR, min_m_n - 1, 1, FLA_TL );
    Ae = ATR; 
  } else {
   uplo = FLA_LOWER_TRIANGULAR;
    FLA_Part_2x2( A, &ATL, &ATR,
                     &ABL, &ABR, 1, min_m_n - 1, FLA_TL );
    Ae = ABL;
  }

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );

  // Bidiag test
  {
    FLA_Obj      norm;
    FLA_Bool     apply_scale;

    FLA_Obj_create( realtype, 1,1, 0,0, &norm );
    FLA_Max_abs_value( A, norm );
    apply_scale = FLA_Obj_gt( norm, FLA_OVERFLOW_SQUARE_THRES ); 
    if ( apply_scale ) FLA_Scal( FLA_SAFE_MIN, A );
    FLA_Bidiag_UT( A, TU, TV );
    if ( apply_scale ) FLA_Bidiag_UT_scale_diagonals( FLA_SAFE_INV_MIN, A ); 
    FLA_Obj_free( &norm );
  }

  // Orthonomal basis U, V. 
  FLA_Obj_create( datatype, m, min_m_n, 0,0, &U ); 
  FLA_Obj_create( datatype, n, min_m_n, 0,0, &V ); 

  FLA_Set( FLA_ZERO, U );
  FLA_Set( FLA_ZERO, V );

  FLA_Bidiag_UT_form_U( A, TU, U );
  FLA_Bidiag_UT_form_V( A, TV, V ); 

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

  // U^H U
  FLA_Obj_create( datatype, min_m_n, min_m_n, 0,0, &DU );
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, U, U, FLA_ZERO, DU );

  // V^H V
  FLA_Obj_create( datatype, min_m_n, min_m_n, 0,0, &DV );
  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, V, V, FLA_ZERO, DV );
  
  // Recover the matrix
  FLA_Obj_create( datatype, min_m_n, min_m_n, 0, 0, &B );
  FLA_Set( FLA_ZERO, B );

  // Set B
  FLA_Obj_create( datatype, min_m_n,     1, 0, 0, &d );  
  FLA_Set_diagonal_vector( A,  d );
  FLA_Set_diagonal_matrix( d, B  );
  FLA_Obj_free( &d );

  if ( min_m_n > 1 ) {
    FLA_Obj_create( datatype, min_m_n -1 , 1, 0, 0, &e );  
    FLA_Set_diagonal_vector( Ae, e );
    if ( uplo == FLA_UPPER_TRIANGULAR ) {
      FLA_Part_2x2( B, &ATL, &ATR,
                    &ABL, &ABR, min_m_n - 1, 1, FLA_TL );
      Be = ATR;
    } else {
      FLA_Part_2x2( B, &ATL, &ATR,
                    &ABL, &ABR, 1, min_m_n - 1, FLA_TL );
      Be = ABL;
    }
    FLA_Set_diagonal_matrix( e, Be );
    FLA_Obj_free( &e );
  }

  // Vb := B (V^H)
  FLA_Obj_create_copy_of( FLA_CONJ_TRANSPOSE, V, &Vb );
  FLA_Trmm_external( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
                     FLA_NONUNIT_DIAG, FLA_ONE, B, Vb );

  // A := U Vb
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, U, Vb, FLA_ZERO, A_recovered );

  residual_A    = FLA_Max_elemwise_diff( A_copy, A_recovered );

  if (1) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - TU - ", TU, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - TV - ", TV, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - B - ", B, "% 6.4e", "------");
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
  FLA_Obj_free( &TU );
  FLA_Obj_free( &TV );

  FLA_Obj_free( &B );

  FLA_Obj_free( &U );
  FLA_Obj_free( &V );
  FLA_Obj_free( &Vb );

  FLA_Obj_free( &DU );
  FLA_Obj_free( &DV );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );


  FLA_Finalize_safe( init_result );     
}
