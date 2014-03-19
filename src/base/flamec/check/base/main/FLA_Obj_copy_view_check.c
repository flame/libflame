/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Obj_copy_view_check( FLA_Obj A, FLA_Obj* B )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( A.base );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( B->base );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

