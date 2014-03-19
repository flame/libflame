/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Svd_check( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_svd_type( jobu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_svd_type( jobv );
  FLA_Check_error_code( e_val );

  // FLA_Svd does not allow FLA_SVD_VECTORS_MIN_OVERWRITE
  // for both jobu and jobv as V cannot be overwritten on A.
  // Use FLA_Svd_ext to allow OVERWRITE options.
  if ( jobu == FLA_SVD_VECTORS_MIN_OVERWRITE ||
       jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
    FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

  // Do not check the jobu and jobv OVERWRITE combination.
  //e_val = FLA_Check_valid_svd_type_combination( jobu, jobv );
  //FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( s );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, s );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( s, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_storage( s );
  FLA_Check_error_code( e_val );

  // When jobu is FLA_SVD_VECTORS_NONE, U may be given without a base object allocated. 
  if ( jobu != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( U ) == FALSE )
  {
    e_val = FLA_Check_identical_object_datatype( A, U );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_length_equals( U, FLA_Obj_length( A ) );
    FLA_Check_error_code( e_val );
  
    // No need to be square.
    //e_val = FLA_Check_square( U );
    //FLA_Check_error_code( e_val );
  }
  
  // When jobv is FLA_SVD_VECTORS_NONE, V may be given without a base object allocated. 
  if ( jobv != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( V ) == FALSE )
  {
    e_val = FLA_Check_identical_object_datatype( A, V );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_length_equals( V, FLA_Obj_width( A ) );
    FLA_Check_error_code( e_val );
    
    // No need to be square.
    //e_val = FLA_Check_square( V );
    //FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

