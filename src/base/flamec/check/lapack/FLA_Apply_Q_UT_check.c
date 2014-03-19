/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_storev( storev );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, W );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    //e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( A ) );
    //FLA_Check_error_code( e_val );
    //e_val = FLA_Check_object_width_equals( T, FLA_Obj_min_dim( A ) );
    //FLA_Check_error_code( e_val );

    if ( storev == FLA_COLUMNWISE )
    {
      e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
      FLA_Check_error_code( e_val );
    }
    else // if ( storev == FLA_ROWWISE )
    {
      e_val = FLA_Check_object_length_equals( B, FLA_Obj_width( A ) );
      FLA_Check_error_code( e_val );
    }
  }
  else // if ( side == FLA_RIGHT )
  {
    //e_val = FLA_Check_object_width_equals( T, FLA_Obj_min_dim( A ) );
    //FLA_Check_error_code( e_val );

    if ( storev == FLA_COLUMNWISE )
    {
      e_val = FLA_Check_object_width_equals( B, FLA_Obj_length( A ) );
      FLA_Check_error_code( e_val );
    }
    else // if ( storev == FLA_ROWWISE )
    {
      e_val = FLA_Check_object_width_equals( B, FLA_Obj_width( A ) );
      FLA_Check_error_code( e_val );
    }
  }

  return FLA_SUCCESS;
}

