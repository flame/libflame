/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Merge_1x2_check( FLA_Obj AL, FLA_Obj AR,   FLA_Obj *A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( AL, AR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_adjacent_objects_1x2( AL, AR );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

