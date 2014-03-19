/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_check( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( G );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_complex_object( G );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( G, A );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( G, FLA_Obj_length( A ) - 1 );
    FLA_Check_error_code( e_val );
  }
  else // if ( side == FLA_RIGHT )
  {
    e_val = FLA_Check_object_length_equals( G, FLA_Obj_width( A ) - 1 );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

