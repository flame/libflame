/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_check( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_blas_trans( transa );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_blas_trans( transb );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( isgn );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( isgn );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_isgn_value( isgn );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_square( B );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_sylv_matrix_dims( A, B, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( scale );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_identical_object_precision( C, scale );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

