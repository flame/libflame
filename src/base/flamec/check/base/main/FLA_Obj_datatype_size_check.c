
#include "FLAME.h"

FLA_Error FLA_Obj_datatype_size_check( FLA_Datatype datatype )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

