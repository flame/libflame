/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Invert( FLA_Conj conj, FLA_Obj x )
{
  FLA_Datatype datatype;
  integer          n_elem;
  integer          inc_x;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING ) 
    FLA_Invert_check( conj, x );

  if ( FLA_Obj_has_zero_dim( x ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( x );

  n_elem   = FLA_Obj_vector_dim( x );
  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bl1_sinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bl1_dinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bl1_cinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  { 
    dcomplex *buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bl1_zinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

