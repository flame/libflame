
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_U_check( FLA_Obj A, FLA_Obj T, FLA_Obj U )
{
  FLA_Error e_val;
  dim_t     m_A, n_A;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, U );
  FLA_Check_error_code( e_val );

  // U is not necessary to be square.
  // e_val = FLA_Check_square( U );
  // FLA_Check_error_code( e_val );

  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width( A );

  // Form U has no problem in overwriting on A which contains house holder vectors
  // on the lower triangular of the diagonal or subdiagonal.
  if ( m_A >= n_A ) 
  {
      e_val = FLA_Check_object_width_equals( T, n_A );
      FLA_Check_error_code( e_val );

      e_val = FLA_Check_object_length_equals( U, m_A );
      FLA_Check_error_code( e_val );
  }
  else
  {
      e_val = FLA_Check_object_width_equals( T, m_A );
      FLA_Check_error_code( e_val );

      e_val = FLA_Check_object_length_equals( U, m_A );
      FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

