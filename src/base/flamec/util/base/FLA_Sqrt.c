/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sqrt( FLA_Obj alpha )
{
  FLA_Datatype datatype;
  integer          r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Sqrt_check( alpha );

  datatype = FLA_Obj_datatype( alpha );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    if ( *buff_alpha <= 0.0F || isnan(*buff_alpha) )
      r_val = FLA_FAILURE;
    else
      *buff_alpha = ( float ) sqrt( *buff_alpha );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    if ( *buff_alpha <= 0.0 || isnan(*buff_alpha) )
      r_val = FLA_FAILURE;
    else
      *buff_alpha = ( double ) sqrt( *buff_alpha );
    
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    if ( buff_alpha->real <= 0.0F || isnan(buff_alpha->real) )
      r_val = FLA_FAILURE;
    else
      buff_alpha->real = ( float ) sqrt( buff_alpha->real );
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    if ( buff_alpha->real <= 0.0 || isnan(buff_alpha->real) )
      r_val = FLA_FAILURE;
    else
      buff_alpha->real = ( double ) sqrt( buff_alpha->real );
    
    break;
  }

  }

  return r_val;
}

