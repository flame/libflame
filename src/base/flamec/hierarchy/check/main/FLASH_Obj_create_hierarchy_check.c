/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Obj_create_hierarchy_check( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* elem_sizes_m, dim_t* elem_sizes_n, FLA_Obj flat_matrix, FLA_Obj* H, unsigned long id, dim_t depth_overall, dim_t* depth_sizes_m, dim_t* depth_sizes_n, dim_t* m_offsets, dim_t* n_offsets )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( elem_sizes_m );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( elem_sizes_n );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( H );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( depth_sizes_m );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( depth_sizes_n );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( m_offsets );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( n_offsets );
  FLA_Check_error_code( e_val );

  // A value of depth < 0 should cause an error.

  // Values of m < 1, n < 1 should cause an error. (or < 0?)

  // First depth entries in depth_sizes_m,_n elem_sizes_m,_n m_,n_offsets should be checked; values < 1 should cause error.

  return FLA_SUCCESS;
}

