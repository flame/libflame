#include "FLAME.h"

#define TESTTYPE FLA_COMPLEX
#define REALTYPE FLA_FLOAT

int main( int argc, char** argv ) {
  FLA_Datatype testtype = TESTTYPE;
  //FLA_Datatype realtype = REALTYPE;
  dim_t        m, n, min_m_n;
  FLA_Obj      A;
  FLA_Obj      a1, b1, rL1, rR1;
  FLA_Obj      a2, b2, rL2, rR2;
  FLA_Uplo     uplo;
  FLA_Error    init_result; 

  if ( argc == 3 ) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    min_m_n = min( m, n );
    uplo = ( m >=n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m uplo\n", argv[0]);
    fprintf(stderr, "       m    : test matrix length\n");
    fprintf(stderr, "       n    : test matrix width\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  // Test matrix A 
  FLA_Obj_create( testtype, m, n, 0, 0, &A );
  FLA_Random_matrix( A );
  {
    scomplex *buff_A = FLA_Obj_buffer_at_view( A );

    buff_A[0].real = -7.9131e-01; buff_A[0].imag =  7.2189e-09; buff_A[2].real = -1.8178e-01; buff_A[2].imag =  1.2624e-01; buff_A[4].real = -5.1535e-02; buff_A[4].imag = -2.2650e-01; buff_A[6].real = -1.0717e-01; buff_A[6].imag = -4.1325e-01;
    buff_A[1].real =  5.9967e-01; buff_A[1].imag = -1.1929e-01; buff_A[3].real =  4.3537e-06; buff_A[3].imag = -2.8598e-06; buff_A[5].real = -2.4200e-01; buff_A[5].imag = -5.6154e-01; buff_A[7].real = -2.7770e-01; buff_A[7].imag =  1.7734e-01;
  }

  FLA_Obj_fshow( stdout,  "- A -", A, "% 6.4e", "--" );

  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &a1 );
  if ( m > 1 ) FLA_Obj_create( testtype, min_m_n-1, 1, 0, 0, &b1 );

  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &rL1 );
  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &rR1 );

  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &a2 );
  if ( m > 1 ) FLA_Obj_create( testtype, min_m_n-1, 1, 0, 0, &b2 );

  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &rL2 );
  FLA_Obj_create( testtype, min_m_n, 1, 0, 0, &rR2 );

  // Mine 
  FLA_Bidiag_UT_extract_diagonals( A, a1, b1 );
  FLA_Obj_fshow( stdout,  "- a1 -", a1, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- b1 -", b1, "% 6.4e", "--" );  

  FLA_Bidiag_UT_realify_diagonals( uplo, a1, b1, rL1, rR1 ); 
  FLA_Obj_fshow( stdout,  "- rL1 -", rL1, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- rR1 -", rR1, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- a1 realified -", a1, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- b1 realified -", b1, "% 6.4e", "--" );  
  
  // Field
  FLA_Bidiag_UT_realify( A, rL2, rR2 );
  FLA_Bidiag_UT_extract_diagonals( A, a2, b2 );
  FLA_Obj_fshow( stdout,  "- rL2 -", rL2, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- rR2 -", rR2, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- a2 realified -", a2, "% 6.4e", "--" );  
  FLA_Obj_fshow( stdout,  "- b2 realified -", b2, "% 6.4e", "--" );  

  printf(" diff_a  = %e\n", FLA_Max_elemwise_diff( a1, a2 ));
  printf(" diff_b  = %e\n", FLA_Max_elemwise_diff( b1, b2 ));
  printf(" diff_rL = %e\n", FLA_Max_elemwise_diff( rL1, rL2 ));
  printf(" diff_rR = %e\n", FLA_Max_elemwise_diff( rR1, rR2 ));

  FLA_Obj_fshow( stdout,  "- A realified-", A, "% 6.4e", "--" );

  FLA_Obj_free( &rR2 );
  FLA_Obj_free( &rL2 );

  FLA_Obj_free( &rR1 );
  FLA_Obj_free( &rL1 );

  if ( m > 1 ) FLA_Obj_free( &b2 );
  if ( m > 1 ) FLA_Obj_free( &b1 );

  FLA_Obj_free( &a2 );
  FLA_Obj_free( &a1 );

  FLA_Obj_free( &A );

  FLA_Finalize_safe( init_result );     
}
