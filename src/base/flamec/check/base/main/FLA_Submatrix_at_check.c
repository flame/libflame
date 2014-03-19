/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Submatrix_at_check( FLA_Datatype datatype, void* buffer, dim_t i, dim_t j, dim_t rs, dim_t cs )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( buffer );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

