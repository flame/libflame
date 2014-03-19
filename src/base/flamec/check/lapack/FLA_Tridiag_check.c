/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_vector( t );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_storage( t );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim_min( t, FLA_Obj_length( A ) - 1 );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

