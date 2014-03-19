/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"
#define EPS 1.e-10

typedef float testtype;
#define TESTTYPE FLA_FLOAT

int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Obj      A, A_copy, p;
  int          *p_buf, i, j, option;
  dim_t        m = 4;
  FLA_Error    init_result; 
  double       residual_A;
  testtype     *A_buf;

  if ( argc == 2 ) 
    option = atoi(argv[1]);
  else 
    option = 1;

  FLA_Init_safe( &init_result );          

  // FLAME Pivot setup
  FLA_Obj_create( datatype, m,m, 0,0, &A );
  FLA_Obj_create( FLA_INT, m,1, 0,0, &p );
  
  // Rand A and create A_copy.
  //FLA_Random_matrix( A ); 
  A_buf = (testtype*)FLA_Obj_buffer_at_view( A );
  for (j=0;j<m;++j)
    for (i=0;i<m;++i)
      A_buf[i+j*m] = (i+1)*10 + (j+1);

  p_buf = (int*)FLA_Obj_buffer_at_view( p );
  p_buf[0] = 2;
  p_buf[1] = 0;
  p_buf[2] = 1;
  p_buf[3] = 0;

  //for (i=0;i<m;++i)
  //  p_buf[i] = (m - i) - 1;
  
  //FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );   
  //FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );   
  //FLA_Obj_fshow( stdout, " - p - ", p, "%d", "------");
    
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_copy );
  FLA_Copy( A, A_copy );

  // Pivot LN + LT
  //option  = 2;

  fprintf(stdout, " Option selected = %d\n", option);
  switch (option) {
  case 1:
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE,  p, A );
    //FLA_Apply_pivots( FLA_LEFT, FLA_TRANSPOSE,     p, A );
    break;
  case 2:
    FLA_Apply_pivots( FLA_RIGHT, FLA_NO_TRANSPOSE, p, A );
    FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE,    p, A );
    break;
  case 3:
    FLA_Apply_pivots( FLA_LEFT,  FLA_NO_TRANSPOSE, p, A );
    FLA_Apply_pivots( FLA_RIGHT, FLA_TRANSPOSE,    p, A );
    FLA_Apply_pivots( FLA_LEFT,  FLA_TRANSPOSE,    p, A );
    FLA_Apply_pivots( FLA_RIGHT, FLA_NO_TRANSPOSE, p, A );
    break;
  }

  // Comapre (A_copy, A_recovered), (y,z) and (y,w)
  residual_A      = FLA_Max_elemwise_diff( A, A_copy );

  if (1 || residual_A > EPS) {
    FLA_Obj_fshow( stdout, " - Given - ", A_copy, "% f", "------");
    FLA_Obj_fshow( stdout, " - Pivoted - ", A, "% f", "------");
    FLA_Obj_fshow( stdout, " - p - ", p, "%d", "------");
    fprintf( stdout, "elem diff A = %12.10e\n\n", residual_A); 
  }
  
  FLA_Obj_free( &A );
  FLA_Obj_free( &A_copy );
  FLA_Obj_free( &p );
  
  FLA_Finalize_safe( init_result );     
}
