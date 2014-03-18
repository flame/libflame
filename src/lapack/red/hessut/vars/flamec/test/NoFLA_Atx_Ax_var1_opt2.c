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

