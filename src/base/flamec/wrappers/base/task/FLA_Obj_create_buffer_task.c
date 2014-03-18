
#include "FLAME.h"

FLA_Error FLA_Obj_create_buffer_task( dim_t rs, dim_t cs, FLA_Obj obj, void* cntl )
{
  FLA_Error r_val;

  r_val = FLA_Obj_create_buffer( rs, cs, &obj );

  FLA_Set( FLA_ZERO, obj );

  return r_val;
}
