/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Repart_2x2_to_3x3_check( FLA_Obj ATL, FLA_Obj ATR,  FLA_Obj *A00, FLA_Obj *A01, FLA_Obj *A02,
                                                                  FLA_Obj *A10, FLA_Obj *A11, FLA_Obj *A12,
                                       FLA_Obj ABL, FLA_Obj ABR,  FLA_Obj *A20, FLA_Obj *A21, FLA_Obj *A22,
                                       dim_t   mb,  dim_t    nb,  FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( ATL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( ABL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( ATR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( ABR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A00 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A10 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A20 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A01 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A02 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  if      ( quadrant == FLA_TL )
  {
    e_val = FLA_Check_attempted_repart_2x2( ATL, mb, nb );
    FLA_Check_error_code( e_val );
  }
  else if ( quadrant == FLA_TR )
  {
    e_val = FLA_Check_attempted_repart_2x2( ATR, mb, nb );
    FLA_Check_error_code( e_val );
  }
  else if ( quadrant == FLA_BL )
  {
    e_val = FLA_Check_attempted_repart_2x2( ABL, mb, nb );
    FLA_Check_error_code( e_val );
  }
  else if ( quadrant == FLA_BR )
  {
    e_val = FLA_Check_attempted_repart_2x2( ABR, mb, nb );
    FLA_Check_error_code( e_val );
  }

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

