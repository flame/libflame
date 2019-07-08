/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Cont_with_3x3_to_2x2_check_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                       FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                          FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                       FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( ATL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ATR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A00 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A10 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A20 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A01 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A02 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Cont_with_3x3_to_2x2_check( FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                       FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                          FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                       FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( ATL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ATR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A00 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A10 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A20 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A01 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A02 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

