
#include "FLAME.h"

FLA_Error FLA_Axpy_buffer_to_object_check( FLA_Trans trans, FLA_Obj alpha, dim_t m, dim_t n, void* A_buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_real_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( B, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A_buffer );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_dims( trans, m, n, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_strides( m, n, rs, cs );
  FLA_Check_error_code( e_val );

  if ( trans == FLA_NO_TRANSPOSE )
  {
    e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, B );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_submatrix_dims_and_offset( n, m, i, j, B );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

