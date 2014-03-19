/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fill_with_linear_dist_check( FLA_Obj shift, FLA_Obj delta, FLA_Obj x )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( shift, delta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( x, delta );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

