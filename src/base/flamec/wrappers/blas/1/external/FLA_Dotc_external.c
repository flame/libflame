/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Dotc_external( FLA_Conj conj, FLA_Obj x, FLA_Obj y, FLA_Obj rho )
{
  FLA_Datatype datatype;
  integer          num_elem;
  integer          inc_x;
  integer          inc_y;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Dotc_check( conj, x, y, rho );

  if ( FLA_Obj_has_zero_dim( x ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  inc_y    = FLA_Obj_vector_inc( y );
  num_elem = FLA_Obj_vector_dim( x );

  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x   = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_y   = ( float * ) FLA_FLOAT_PTR( y );
    float *buff_rho = ( float * ) FLA_FLOAT_PTR( rho );

    bl1_sdot( blis_conj,
              num_elem, 
              buff_x, inc_x,
              buff_y, inc_y,
              buff_rho );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x   = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_y   = ( double * ) FLA_DOUBLE_PTR( y );
    double *buff_rho = ( double * ) FLA_DOUBLE_PTR( rho );

    bl1_ddot( blis_conj,
              num_elem, 
              buff_x, inc_x,
              buff_y, inc_y,
              buff_rho );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x   = ( scomplex * ) FLA_COMPLEX_PTR( x );
    scomplex *buff_y   = ( scomplex * ) FLA_COMPLEX_PTR( y );
    scomplex *buff_rho = ( scomplex * ) FLA_COMPLEX_PTR( rho );

    bl1_cdot( blis_conj,
              num_elem, 
              buff_x, inc_x,
              buff_y, inc_y,
              buff_rho );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_x   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    dcomplex *buff_y   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( y );
    dcomplex *buff_rho = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( rho );

    bl1_zdot( blis_conj,
              num_elem, 
              buff_x, inc_x,
              buff_y, inc_y,
              buff_rho );

    break;
  }

  }

  return FLA_SUCCESS;
}

