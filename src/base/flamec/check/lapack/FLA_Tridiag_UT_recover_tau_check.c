/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_recover_tau_check( FLA_Obj T, FLA_Obj tau )
{
  FLA_Error    e_val;

  e_val = FLA_Check_floating_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( T, tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( T, FLA_Obj_vector_dim( tau ) + 1 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

