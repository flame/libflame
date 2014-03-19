/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Repart_1x2_to_1x3_check( FLA_Obj  AL,              FLA_Obj  AR,
                                       FLA_Obj *A0, FLA_Obj *A1, FLA_Obj *A2,
                                                    dim_t    nb, FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A0 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  if      ( side == FLA_LEFT )
  {
    e_val = FLA_Check_attempted_repart_1x2( AL, nb );
    FLA_Check_error_code( e_val );
  }
  else if ( side == FLA_RIGHT )
  {
    e_val = FLA_Check_attempted_repart_1x2( AR, nb );
    FLA_Check_error_code( e_val );
  }

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

