/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Mach_params_check( FLA_Machval machval, FLA_Obj val )
{
  FLA_Error    e_val;

  e_val = FLA_Check_valid_machval( machval );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( val );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

