/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR2_UT_check( FLA_Obj B, FLA_Obj D, FLA_Obj T )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( B, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( B, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, B, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, T, B, T );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

