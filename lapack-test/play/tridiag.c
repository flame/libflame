#include "FLAME.h"

//typedef float testtype;
//#define TESTTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE

typedef dcomplex testtype;
#define TESTTYPE FLA_DOUBLE_COMPLEX

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Obj      A, T, A_copy, A_recovered;     
  FLA_Obj      Q, Qb, B;
  FLA_Uplo     uplo;
  dim_t        m;
  FLA_Error    init_result; 
  double       residual = 0.0;

  if ( argc == 3 ) {
    m = atoi( argv[1] );
    if ( *argv[2] == 'u' || *argv[2] == 'U' )
      uplo = FLA_UPPER_TRIANGULAR;
    else 
      uplo = FLA_LOWER_TRIANGULAR;
  } else if ( argc == 2 ) {
    m = atoi( argv[1] );
    uplo = FLA_LOWER_TRIANGULAR;
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m [uplo] \n", argv[0]);
    fprintf(stderr, "       m    : matrix length (square)\n");
    fprintf(stderr, "       uplo : u || U (upper), default is lower \n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 ) 
    return 0;

  FLA_Init_safe( &init_result );          

  // FLAME Tridiag setup
  FLA_Obj_create( datatype, m, m, 0, 0, &A );
  FLA_Tridiag_UT_create_T( A, &T );

  // Rand A and create A_copy.
  FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A ); 
  FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_recovered );

  // Tridiag test
  FLA_Tridiag_UT( uplo, A, T );

  // Orthonomal basis Q (lower) or Q^H (upper) 
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &Q );  
  FLA_Tridiag_UT_form_Q ( uplo, Q, T, Q );

  if ( FLA_Obj_is_complex( A ) ){
    FLA_Obj r;

    FLA_Obj_create( datatype, m, 1, 0, 0, &r );
    FLA_Tridiag_UT_realify( uplo, A, r );

    FLA_Obj_fshow( stdout, " - r - ", r, "% 6.4e", "------");

    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE,    r, Q );

    FLA_Obj_free( &r );
  }

  // Q'Q
  {
    FLA_Obj DQ, norm;
    FLA_Obj_create( datatype, m, m, 0, 0, &DQ );
    FLA_Set_to_identity( DQ );
    FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, Q, Q, FLA_MINUS_ONE, DQ );

    FLA_Obj_create( datatype, 1, 1, 0, 0, &norm );
    FLA_Max_abs_value( DQ, norm );
    FLA_Obj_fshow( stdout, " - Q^H Q - I  - ", norm, "% 6.4e", "------");
    FLA_Obj_free( &norm );
    FLA_Obj_free ( &DQ );
  }

  // Recover the matrix
  FLA_Obj_create( datatype, m, m, 0, 0, &B );
  FLA_Set( FLA_ZERO, B );

  // Set B
  {
    FLA_Obj ATL, ATR,
            ABL, ABR, Ae;
    FLA_Obj d, e;

    if ( uplo == FLA_UPPER_TRIANGULAR ) {
      FLA_Part_2x2( A, &ATL, &ATR,
                       &ABL, &ABR, 1, 1, FLA_BL );
      Ae = ATR;
    } else {
      FLA_Part_2x2( A, &ATL, &ATR,
                       &ABL, &ABR, 1, 1, FLA_TR );
      Ae = ABL;
    }

    FLA_Obj_create( datatype, m, 1, 0, 0, &d );  
    FLA_Set_diagonal_vector( A, d );
    FLA_Set_diagonal_matrix( d, B );
    FLA_Obj_free( &d );

    if ( m > 1 ) {
      FLA_Obj_create( datatype, m -1 , 1, 0, 0, &e );    
      FLA_Set_diagonal_vector( Ae, e );
      Ae.base = B.base;
      FLA_Set_diagonal_matrix( e, Ae );
      FLA_Obj_free( &e );

      FLA_Hermitianize( uplo, B );
    }
  }

  // Recover A from Q B Q^H
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, Q, &Qb );

  // Qb := B (Q^H)
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                     FLA_ONE, B, Q, FLA_ZERO, Qb );
  // A := Q Qb
  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_ONE, Q, Qb, FLA_ZERO, A_recovered );
  
  residual = FLA_Max_elemwise_diff( A_copy, A_recovered );

  if (1) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Factor - ", A, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - T - ", T, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - B - ", B, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Q - ", Q, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Qb - ", Qb, "% 6.4e", "------");
    FLA_Obj_fshow( stdout, " - Recovered A - ", A_recovered, "% 6.4e", "------");
    fprintf( stdout, "lapack2flame: %lu x %lu: ", m, m);
    fprintf( stdout, "recovery A = %12.10e\n\n", residual );
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &T );

  FLA_Obj_free( &B );

  FLA_Obj_free( &Q );
  FLA_Obj_free( &Qb );

  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &A_recovered );

  FLA_Finalize_safe( init_result );     
}
