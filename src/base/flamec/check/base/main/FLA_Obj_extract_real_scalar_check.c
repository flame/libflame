/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Obj_extract_real_scalar_check( FLA_Obj alpha, double* alpha_value )
{
  FLA_Error e_val;

  e_val = FLA_Check_real_object( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( alpha_value );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}


