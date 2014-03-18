
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_extract_diagonals_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Error e_val;
  dim_t     m_A;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  m_A   = FLA_Obj_length( A );

  e_val = FLA_Check_nonconstant_object( d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( d, m_A );
  FLA_Check_error_code( e_val );

  if ( m_A > 1 ) 
  {
    e_val = FLA_Check_nonconstant_object( e );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_identical_object_datatype( A, e );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_if_vector( e );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_vector_dim( e, m_A - 1 );
    FLA_Check_error_code( e_val );
  }
  
  return FLA_SUCCESS;
}

