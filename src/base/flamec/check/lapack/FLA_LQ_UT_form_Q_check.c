
#include "FLAME.h"

FLA_Error FLA_LQ_UT_form_Q_check( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, Q );
  FLA_Check_error_code( e_val );

  // - the width of T represents the number of applied house-holder vectors.
  // - the length of A may represent the same number of those vectors.
  // - however, if A and Q share the same buffer (in-place operation),
  //   then the length of A should match to the length of the maximum length 
  //   of the applied house-holder vectors.
  //if ( FLA_Obj_is_overlapped( A, Q ) == FALSE ) 
  //{
  //  e_val = FLA_Check_object_width_equals( T, FLA_Obj_length( A ) );
  //  FLA_Check_error_code( e_val );
  //}

  e_val = FLA_Check_object_width_equals( Q, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  // Q matrix should not be restricted to be a square matrix
  // e_val = FLA_Check_square( Q );
  // FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

