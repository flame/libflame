/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_UDdate_UT_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( C ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( D ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

