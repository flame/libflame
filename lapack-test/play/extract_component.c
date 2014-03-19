/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#define TESTTYPE FLA_DOUBLE
#define REALTYPE FLA_DOUBLE

int main( int argc, char** argv ) {
  FLA_Datatype testtype = TESTTYPE;
  FLA_Datatype realtype = REALTYPE;
  dim_t        m;
  FLA_Obj      a, b;
  FLA_Error    init_result; 

  if ( argc == 2 ) {
    m = atoi(argv[1]);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m\n", argv[0]);
    fprintf(stderr, "       m       : test vector length\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          
  
  FLA_Obj_create( testtype, m, 1, 0, 0, &a );
  FLA_Random_matrix( a );
  FLA_Obj_fshow( stdout,  "- a -", a, "% 6.4e", "--" );

  FLA_Obj_create( realtype, 1, m, 0, 0, &b );
  
  FLA_Obj_extract_real_part( a, b );
  FLA_Obj_fshow( stdout,  "- a real -", b, "% 6.4e", "--" );

  FLA_Obj_extract_imag_part( a, b );
  FLA_Obj_fshow( stdout,  "- a imag -", b, "% 6.4e", "--" );

  FLA_Obj_free( &b );
  FLA_Obj_free( &a );

  FLA_Finalize_safe( init_result );     
}
