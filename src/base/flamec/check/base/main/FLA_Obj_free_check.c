/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
/* Not Required */

FLA_Error FLA_Obj_free_check( FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  // FLA_Free check base pointer; thus, this is not necessary.
  //e_val = FLA_Check_null_pointer( obj->base );
  //FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

