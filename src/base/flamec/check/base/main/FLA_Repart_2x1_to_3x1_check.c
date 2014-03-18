
#include "FLAME.h"

FLA_Error FLA_Repart_2x1_to_3x1_check( FLA_Obj AT,  FLA_Obj *A0,
                                                    FLA_Obj *A1,
                                       FLA_Obj AB,  FLA_Obj *A2,
                                       dim_t   mb,  FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A0 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_topbottom_side( side );
  FLA_Check_error_code( e_val );

  if      ( side == FLA_TOP )
  {
    e_val = FLA_Check_attempted_repart_2x1( AT, mb );
    FLA_Check_error_code( e_val );
  }
  else if ( side == FLA_BOTTOM )
  {
    e_val = FLA_Check_attempted_repart_2x1( AB, mb );
    FLA_Check_error_code( e_val );
  }

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

