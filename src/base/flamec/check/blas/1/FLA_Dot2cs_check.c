/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Dot2cs_check( FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj beta, FLA_Obj rho )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_conj( conj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( x, y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( x, rho );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( x, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( x, beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( rho );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_equal_vector_dims( x, y );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

