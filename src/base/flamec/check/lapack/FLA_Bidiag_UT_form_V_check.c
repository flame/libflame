
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_form_V_check( FLA_Obj A, FLA_Obj T, FLA_Obj V )
{
  FLA_Error e_val;
  dim_t     m_A, n_A;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, V );
  FLA_Check_error_code( e_val );

  // U is not necessary to be square.
  // e_val = FLA_Check_square( V );
  // FLA_Check_error_code( e_val );
  

  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width( A );

  // Form V (not V^H) has a problem that dimensions are mismatched 
  // when it overwrites on A that contains house holder vectors on
  // the uppper triangular of the diagonal or subdiagonal.
  if ( m_A >= n_A )
  {
    e_val = FLA_Check_object_width_equals( T, n_A );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_object_width_equals( V, n_A );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_object_width_equals( T, m_A );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_object_length_equals( V, n_A );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

