/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Svd_compute_scaling_check( FLA_Obj A, FLA_Obj sigma )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, sigma );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

