/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_FS_incpiv_check( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, L );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_nonconstant_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, L );
  FLA_Check_error_code( e_val );

  //e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, A, B, B );
  //FLA_Check_error_code( e_val );

  // Until we update FS_incpiv to support multiple right-hand sides, we force B
  // to have only one column (ie: we force B to be a vector).
  e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( B, 1 );
  FLA_Check_error_code( e_val );


  return FLA_SUCCESS;
}

