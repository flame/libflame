/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Obj_create_ext_check( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t m, dim_t n, dim_t m_inner, dim_t n_inner, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_elemtype( elemtype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  // If both m and n are zero, we do not need to check cs/rs.
  if ( m > 0 && n > 0 )
  {
    e_val = FLA_Check_matrix_strides( m, n, rs, cs );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

