/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Obj_equals_check_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype_ts( FLA_cntl_init_i, A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Obj_equals_check( FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

