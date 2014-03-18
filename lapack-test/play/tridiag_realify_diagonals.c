#include "FLAME.h"

#define TESTTYPE FLA_COMPLEX

int main( int argc, char** argv ) {
  FLA_Datatype testtype = TESTTYPE;
  dim_t        m;
  FLA_Obj      A;
  FLA_Obj      a1, b1, r1;
  FLA_Obj      a2, b2, r2;
  FLA_Uplo     uplo;
  FLA_Error    init_result; 

  if ( argc == 3 ) {
    m = atoi(argv[1]);
    uplo = ( atoi(argv[2]) == 1  ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m uplo\n", argv[0]);
    fprintf(stderr, "       m    : test matrix length\n");
    fprintf(stderr, "       uplo : 0) lower, 1) upper\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // Test matrix A 
  FLA_Obj_create( testtype, m, m, 0, 0, &A );
  FLA_Random_spd_matrix( uplo, A );
  FLA_Hermitianize( uplo, A );
  FLA_Obj_fshow( stdout,  "- A -", A, "% 6.4e", "--" );

  FLA_Obj_create( testtype, m, 1, 0, 0, &a1 );
  FLA_Obj_create( testtype, m, 1, 0, 0, &a2 );

  if ( m > 1 ) {
    FLA_Obj_create( testtype, m-1, 1, 0, 0, &b1 );
    FLA_Obj_create( testtype, m-1, 1, 0, 0, &b2 );
  }
  
  FLA_Obj_create( testtype, m, 1, 0, 0, &r1 );
  FLA_Obj_create( testtype, m, 1, 0, 0, &r2 );

  // Mine 
  FLA_Tridiag_UT_extract_diagonals( uplo, A, a1, b1 );
  FLA_Obj_fshow( stdout,  "- a1 -", a1, "% 6.4e", "--" );  
  if ( m > 1 ) FLA_Obj_fshow( stdout,  "- b1 -", b1, "% 6.4e", "--" );  

  FLA_Tridiag_UT_realify_subdiagonal( b1, r1 );
  if ( m > 1 ) FLA_Obj_fshow( stdout,  "- b1 realified -", b1, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- r1 -", r1, "% 6.4e", "--" );  

  
  // Field
  FLA_Tridiag_UT_realify( uplo, A, r2 );
  FLA_Tridiag_UT_extract_diagonals( uplo, A, a2, b2 );
  FLA_Obj_fshow( stdout,  "- a2  -", a2, "% 6.4e", "--" );  
  if ( m > 1 ) FLA_Obj_fshow( stdout,  "- b2 realified -", b2, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- r2 -", r2, "% 6.4e", "--" );  

  printf(" diff_a  = %e\n", FLA_Max_elemwise_diff( a1, a2 ));
  if ( m > 1 ) printf(" diff_b  = %e\n", FLA_Max_elemwise_diff( b1, b2 ));
  printf(" diff_rL = %e\n", FLA_Max_elemwise_diff( r1, r2 ));

  FLA_Obj_fshow( stdout,  "- A realified-", A, "% 6.4e", "--" );

  FLA_Obj_free( &r2 );
  FLA_Obj_free( &r1 );

  if ( m > 1 ) {
    FLA_Obj_free( &b2 );
    FLA_Obj_free( &b1 );
  }

  FLA_Obj_free( &a2 );
  FLA_Obj_free( &a1 );

  FLA_Obj_free( &A );

  FLA_Finalize_safe( init_result );     
}
