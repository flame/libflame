/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_diag_matrix_check( FLA_Side side, FLA_Conj conj, FLA_Obj x, FLA_Obj A )
{
  FLA_Error    e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_conj( conj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, x );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( A, FLA_Obj_vector_dim( x ) );
    FLA_Check_error_code( e_val );
  }
  else // if ( side == FLA_RIGHT )
  {
    e_val = FLA_Check_object_width_equals( A, FLA_Obj_vector_dim( x ) );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

