
#include "FLAME.h"

FLA_Error FLA_Obj_free_buffer_task( FLA_Obj obj, void* cntl )
{
  return FLA_Obj_free_buffer( &obj );
}
