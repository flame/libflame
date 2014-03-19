/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Form_perm_matrix_check( FLA_Obj p, FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  FLA_Check_matrix_vector_dims( FLA_NO_TRANSPOSE, A, p, p );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

