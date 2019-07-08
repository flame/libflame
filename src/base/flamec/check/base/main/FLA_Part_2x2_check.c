/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Part_2x2_check_ts( FLA_cntl_init_s *FLA_cntl_init_i,
                                 FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                          FLA_Obj *A21, FLA_Obj *A22,
                              dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype_ts( FLA_cntl_init_i, A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Part_2x2_check( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                          FLA_Obj *A21, FLA_Obj *A22,
                              dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

