/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Obj_fshow_check( FILE* file, char* s1, FLA_Obj obj, char* format, char* s2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( file );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( s1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_scalar_elemtype( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( format );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( s2 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

