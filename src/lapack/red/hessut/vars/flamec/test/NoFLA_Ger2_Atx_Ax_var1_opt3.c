/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
  //// MyFLA_Obj_set_to_zero( w );
  for ( i = 0; i < m; i++ ) {
    ww( i ) = 0.0;
  }

  //// m = FLA_Obj_length( A );
  //// n = FLA_Obj_width ( A );
  //// ldA = FLA_Obj_ldim( A );
  
  //// buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  //// buff_u = ( double * ) FLA_Obj_buffer_at_view( u );
  //// buff_y = ( double * ) FLA_Obj_buffer_at_view( y );
  //// buff_z = ( double * ) FLA_Obj_buffer_at_view( z );
  //// buff_x = ( double * ) FLA_Obj_buffer_at_view( x );
  //// buff_v = ( double * ) FLA_Obj_buffer_at_view( v );
  //// buff_w = ( double * ) FLA_Obj_buffer_at_view( w );

  //// printf( " pepe\n" );
 
  for ( j = 0; j < n; j++ ) {
    temp1 = *beta * yy( j );
    temp2 = *beta * buff_u[ j ];
    // daxpy_( &m, &temp1, buff_u, &i_one, &AA( 0,j ), &i_one );
    // daxpy_( &m, &temp2, buff_z, &i_one, &AA( 0, j ), &i_one );
    buff_A_col_j = &( buff_A[ 0 + ldA * j ] );
    pA = buff_A_col_j;
    pu = buff_u;
    pz = buff_z;
    for( i = 0; i < m; i++ ) {
      temp3 = temp1 * *pu;
      *pA += temp3 + temp2 * *pz;
      pA++;
      pu++;
      pz++;
    }
    //// daxpy_( &m, &uu( j ), buff_z, &i_one, &AA( 0,j ), &i_one );

    vv( j ) =  ddot_( &m, &AA( 0, j ), &i_one, buff_x, &i_one );
    //bl1_ddot( BLIS1_CONJUGATE,
    //          m,
    //          &AA(0,j), i_one,
    //          buff_x,   i_one,
    //          &vv(j) );
    daxpy_( &m, &xx( j ), &AA( 0, j ), &i_one, buff_w, &i_one );
    //bl1_daxpyv( BLIS1_NO_CONJUGATE,
    //            m,
    //            &xx(j),
    //            &AA(0,j), i_one,
    //            buff_w,   i_one );

  }  

  return FLA_SUCCESS;
}

