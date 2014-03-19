/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Obj_attach_buffer_check( void *buffer, dim_t rs, dim_t cs, FLA_Obj* H )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( H );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_strides( FLASH_Obj_base_scalar_length( *H ), FLASH_Obj_base_scalar_width( *H ), rs, cs );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

