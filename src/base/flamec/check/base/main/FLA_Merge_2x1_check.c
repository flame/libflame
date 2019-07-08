/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Merge_2x1_check_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj AT,
                               FLA_Obj AB,  FLA_Obj *A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( AT, AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_adjacent_objects_2x1( AT,
                                          AB );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Merge_2x1_check( FLA_Obj AT,
                               FLA_Obj AB,  FLA_Obj *A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( AT, AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_adjacent_objects_2x1( AT,
                                          AB );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

