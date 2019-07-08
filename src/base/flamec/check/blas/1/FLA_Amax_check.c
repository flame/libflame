/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Amax_check_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj x, FLA_Obj index )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object_ts( FLA_cntl_init_i, x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object_ts( FLA_cntl_init_i, x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( index ); 
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object_ts( FLA_cntl_init_i, index );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( index );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Amax_check( FLA_Obj x, FLA_Obj index )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( index ); 
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( index );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( index );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

