/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
  //// MyFLA_Obj_set_to_zero( v );
  //// MyFLA_Obj_set_to_zero( w );
  for ( j = 0; j < n; j++ ) {
    vv( j ) = 0.0;
    ww( j ) = 0.0;
  }

  for ( j = 0; j < n; j++ ) {
    vv( j ) =  ddot_( &m, &AA( 0,j ), &i_one, buff_x, &i_one );
    daxpy_( &m, &xx( j ), &AA( 0,j ), &i_one, buff_w, &i_one );
  }  

  return FLA_SUCCESS;
}

