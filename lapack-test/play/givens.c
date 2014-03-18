#include "FLAME.h"

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE
//#define COMPTYPE FLA_DOUBLE_COMPLEX

typedef float testtype;
#define TESTTYPE FLA_FLOAT
#define REALTYPE FLA_FLOAT
#define COMPTYPE FLA_COMPLEX

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Datatype comptype = COMPTYPE;
  dim_t        m;
  FLA_Obj      eye, A, G, A_copy, norm;
  FLA_Obj      GT, G0,
                   G1,
               GB, G2;
  FLA_Error    init_result; 
  double       residual;

  if ( argc == 2 ) {
    m = atoi(argv[1]);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m flip \n", argv[0]);
    fprintf(stderr, "       m    : matrix length\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          
  
  FLA_Obj_create( comptype, m-1, 1, 0, 0, &G );
  FLA_Random_matrix( G );

  FLA_Obj_create( comptype, 1, 1, 0, 0, &norm );

  // Normalization
  FLA_Part_2x1( G,    &GT,  
                      &GB,      0, FLA_TOP );
  while ( FLA_Obj_length( GB ) > 0 ) {
    FLA_Repart_2x1_to_3x1( GT, &G0,
                               &G1,
                           GB, &G2, 1, FLA_BOTTOM );
    // --------------------------------------------
    FLA_Copy( G1, norm );
    FLA_Absolute_value( norm );
    FLA_Inv_scal( norm, G1 );
    // --------------------------------------------
    FLA_Cont_with_3x1_to_2x1( &GT, G0,
                                   G1,
                              &GB, G2, FLA_TOP );
  }
  FLA_Obj_free( &norm );

  FLA_Obj_fshow( stdout, " - Rand G - ", G, "% 6.4e", "--");
  
  FLA_Obj_create( datatype, m, m, 0, 0, &A );
  FLA_Random_matrix( A );

  //FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
  //FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );



  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_copy );

  {
    fprintf( stdout, " == Row-major and column-major Givens test ==\n");

    // right forward, 
    FLA_Obj_fshow( stdout, " - A col - ", A, "% 6.4e", "--");
    FLA_Obj_flip_base( &A );
    FLA_Copy( A_copy, A );
    FLA_Obj_fshow( stdout, " - A row - ", A, "% 6.4e", "--");
    
    FLA_Apply_G_rf_blk_var3( G, A,      m ); // row-major
    FLA_Apply_G_rf_blk_var3( G, A_copy, m ); // column-major

    FLA_Obj_fshow( stdout, " - AG row - ", A,      "% 6.4e", "--");
    FLA_Obj_fshow( stdout, " - AG col - ", A_copy, "% 6.4e", "--");

    // this should be same
    residual = FLA_Max_elemwise_diff( A_copy, A );
    fprintf( stdout, "lapack2flame: %lu x %lu: ", m, m );
    fprintf( stdout, "| A - A_copy | = %12.10e \n\n", residual ) ;

    // A is row major afterwards, A_copy is column major.
  }

  FLA_Obj_create( datatype, m, m, 0, 0, &eye );
  FLA_Set_to_identity( A );
  FLA_Set_to_identity( A_copy );
  FLA_Set_to_identity( eye );
  {
    fprintf( stdout, " == Row-major and column-major Givens test ==\n");
    FLA_Apply_G_rf_blk_var3( G, A, m ); // eye G
    FLA_Apply_G_lf_blk_var3( G, A, m ); // G^H eye G = eye

    FLA_Apply_G_rf_blk_var3( G, A_copy, m ); // eye G
    FLA_Apply_G_lf_blk_var3( G, A_copy, m ); // G^H eye G = eye
    
    // this should be eye
    fprintf( stdout, "lapack2flame: %lu x %lu\n", m, m );
    residual = FLA_Max_elemwise_diff( A, eye );
    fprintf( stdout, "(row major)    | G^H G - eye | = %12.10e \n\n", residual ) ;
    residual = FLA_Max_elemwise_diff( A_copy, eye );
    fprintf( stdout, "(column major) | G^H G - eye | = %12.10e \n\n", residual ) ;
  }
  FLA_Obj_free( &eye );    
  FLA_Obj_free( &G );
  FLA_Obj_free( &A );
  FLA_Obj_free( &A_copy );

  FLA_Finalize_safe( init_result );     
}
