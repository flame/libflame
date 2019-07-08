/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_LU_piv_check_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj A, FLA_Obj p )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object_ts( FLA_cntl_init_i, A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object_ts( FLA_cntl_init_i, A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object_ts( FLA_cntl_init_i, p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_vector( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim_min( p, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_LU_piv_check( FLA_Obj A, FLA_Obj p )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_vector( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim_min( p, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

