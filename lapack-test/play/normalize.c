#include "FLAME.h"

typedef float testtype;
#define TESTTYPE FLA_FLOAT
#define REALTYPE FLA_FLOAT
#define COMPTYPE FLA_COMPLEX

int main( int argc, char** argv ) {
  FLA_Datatype comptype = COMPTYPE;
  FLA_Datatype realtype = REALTYPE;
  dim_t        m;
  FLA_Obj      a, aT, aB, a0, a1, a2;
  FLA_Obj      v, vT, vB, v0, v1, v2;
  FLA_Error    init_result; 
  int          use_abs = 1;

  if ( argc == 3 ) {
    m = atoi(argv[1]);
    use_abs = atoi(argv[2]);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m use_abs\n", argv[0]);
    fprintf(stderr, "       m       : test vector length\n");
    fprintf(stderr, "       use_abs : 0 - norm (realtype), 1 - abs (complex type)\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          
  
  FLA_Obj_create( comptype, m, 1, 0, 0, &a );
  FLA_Obj_create( use_abs ? comptype : realtype, m, 1, 0, 0, &v );

  FLA_Random_matrix( a );
  FLA_Set( FLA_ZERO, v );

  FLA_Obj_fshow( stdout,  "- a -", a, "% 6.4e", "--" );

  // Normalize a vector
  FLA_Part_2x1( a,    &aT,  
                      &aB,      0, FLA_TOP );
  FLA_Part_2x1( v,    &vT,  
                      &vB,      0, FLA_TOP );
  while ( FLA_Obj_length( aB ) > 0 ) {
    FLA_Repart_2x1_to_3x1( aT, &a0,
                               &a1,
                           aB, &a2, 1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( vT, &v0,
                               &v1,
                           vB, &v2, 1, FLA_BOTTOM );
    // --------------------------------------------
    if ( use_abs ) { // a and v are complex datatype
      FLA_Copy( a1, v1 );
      FLA_Absolute_value( v1 );
    } else {         // v is real datatype 
      FLA_Nrm2( a1, v1 );
    }
    if ( FLA_Obj_equals( v1, FLA_ZERO ) )
      printf( " ZERO DETECTED\n" );
    else
      FLA_Inv_scal( v1, a1 ); // Normalize the scalar
    // --------------------------------------------
    FLA_Cont_with_3x1_to_2x1( &aT, a0,
                                   a1,
                              &aB, a2, FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &vT, v0,
                                   v1,
                              &vB, v2, FLA_TOP );
  }

  FLA_Obj_fshow( stdout,  "- a -", a, "% 6.4e", "--" );
  FLA_Obj_fshow( stdout,  "- v -", v, "% 6.4e", "--" );

  // Check whether it is normalized
  FLA_Part_2x1( a,    &aT,  
                      &aB,      0, FLA_TOP );
  FLA_Part_2x1( v,    &vT,  
                      &vB,      0, FLA_TOP );
  while ( FLA_Obj_length( aB ) > 0 ) {
    FLA_Repart_2x1_to_3x1( aT, &a0,
                               &a1,
                           aB, &a2, 1, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( vT, &v0,
                               &v1,
                           vB, &v2, 1, FLA_BOTTOM );
    // --------------------------------------------
    if ( use_abs ) { // a and v are same datatype
      FLA_Copy( a1, v1 );
      FLA_Absolute_value( v1 );
    } else {         // v is realdatatype 
      FLA_Nrm2( a1, v1 );
    }
    // --------------------------------------------
    FLA_Cont_with_3x1_to_2x1( &aT, a0,
                                   a1,
                              &aB, a2, FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &vT, v0,
                                   v1,
                              &vB, v2, FLA_TOP );
  }

  FLA_Obj_fshow( stdout, " - all should be one - ", v, "% 6.4e", "--");

  FLA_Obj_free( &a );
  FLA_Obj_free( &v );

  FLA_Finalize_safe( init_result );     
}
