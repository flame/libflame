
#include "FLAME.h"

FLA_Error FLA_Shift_pivots_to_check( FLA_Pivot_type ptype, FLA_Obj p )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_pivot_type( ptype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_vector( p );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

