/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Amax_external( FLA_Obj x, FLA_Obj index )
{
  FLA_Datatype datatype;
  int          num_elem;
  int          inc_x;
  int         *buff_index;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Amax_check( x, index );

  buff_index = ( int * ) FLA_INT_PTR( index );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    *buff_index = 0;
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  num_elem = FLA_Obj_vector_dim( x );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bl1_samax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }
  
  case FLA_DOUBLE:
  {
    double* buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bl1_damax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }
  
  case FLA_COMPLEX:
  {
    scomplex* buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bl1_camax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bl1_zamax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }

  }

  return FLA_SUCCESS;
}

