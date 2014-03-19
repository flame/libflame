/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_form_Q_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
  FLA_Error e_val;
  dim_t     n_T;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  n_T   = FLA_Obj_width( T );
  
  e_val = FLA_Check_object_width_equals( A, n_T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( Q );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( Q );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( Q );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, Q );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( Q, n_T );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

