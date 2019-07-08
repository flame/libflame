/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Cont_with_3x1_to_2x1_check_ts( FLA_cntl_init_s *FLA_cntl_init_i,
                                             FLA_Obj *AT,  FLA_Obj A0,
                                                        FLA_Obj A1,
                                             FLA_Obj *AB,  FLA_Obj A2,
                                                        FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A0 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_topbottom_side( side );
  FLA_Check_error_code( e_val );

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Cont_with_3x1_to_2x1_check( FLA_Obj *AT,  FLA_Obj A0,
                                                        FLA_Obj A1,
                                          FLA_Obj *AB,  FLA_Obj A2,
                                                        FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A0 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_topbottom_side( side );
  FLA_Check_error_code( e_val );

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

