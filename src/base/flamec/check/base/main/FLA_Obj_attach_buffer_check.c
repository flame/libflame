
#include "FLAME.h"

FLA_Error FLA_Obj_attach_buffer_check( void *buffer, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_strides( FLA_Obj_length( *obj ), FLA_Obj_width( *obj ), rs, cs );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

