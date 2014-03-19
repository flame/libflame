/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Wilkshift_tridiag_check( FLA_Obj delta1, FLA_Obj epsilon, FLA_Obj delta2, FLA_Obj kappa )
{
  FLA_Error e_val;

  e_val = FLA_Check_nonconstant_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, epsilon );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, kappa );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( kappa );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

