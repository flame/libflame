/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Dots_external( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj beta, FLA_Obj rho )
{
  FLA_Datatype datatype;
  int          num_elem;
  int          inc_x;
  int          inc_y;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Dots_check( alpha, x, y, beta, rho );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    FLA_Scal_external( beta, rho );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  inc_y    = FLA_Obj_vector_inc( y );
  num_elem = FLA_Obj_vector_dim( x );

  FLA_Param_map_flame_to_blis_conj( FLA_NO_CONJUGATE, &blis_conj );

  switch ( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_x      = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_y      = ( float * ) FLA_FLOAT_PTR( y );
    float *buff_rho    = ( float * ) FLA_FLOAT_PTR( rho );
    float *buff_alpha  = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta   = ( float * ) FLA_FLOAT_PTR( beta );

    bl1_sdots( blis_conj,
               num_elem, 
               buff_alpha,
               buff_x, inc_x, 
               buff_y, inc_y,
               buff_beta,
               buff_rho );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x      = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_y      = ( double * ) FLA_DOUBLE_PTR( y );
    double *buff_rho    = ( double * ) FLA_DOUBLE_PTR( rho );
    double *buff_alpha  = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta   = ( double * ) FLA_DOUBLE_PTR( beta );

    bl1_ddots( blis_conj,
               num_elem, 
               buff_alpha,
               buff_x, inc_x, 
               buff_y, inc_y,
               buff_beta,
               buff_rho );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x      = ( scomplex * ) FLA_COMPLEX_PTR( x );
    scomplex *buff_y      = ( scomplex * ) FLA_COMPLEX_PTR( y );
    scomplex *buff_rho    = ( scomplex * ) FLA_COMPLEX_PTR( rho );
    scomplex *buff_alpha  = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta   = ( scomplex * ) FLA_COMPLEX_PTR( beta );

    bl1_cdots( blis_conj,
               num_elem, 
               buff_alpha,
               buff_x, inc_x, 
               buff_y, inc_y,
               buff_beta,
               buff_rho );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_x      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    dcomplex *buff_y      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( y );
    dcomplex *buff_rho    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( rho );
    dcomplex *buff_alpha  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bl1_zdots( blis_conj,
               num_elem, 
               buff_alpha,
               buff_x, inc_x, 
               buff_y, inc_y,
               buff_beta,
               buff_rho );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

