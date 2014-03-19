/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Obj_create_buffer_task( dim_t rs, dim_t cs, FLA_Obj obj, void* cntl )
{
  FLA_Error r_val;

  r_val = FLA_Obj_create_buffer( rs, cs, &obj );

  FLA_Set( FLA_ZERO, obj );

  return r_val;
}
