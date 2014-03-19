/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_realify_diagonals_check( FLA_Uplo uplo, FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e )
{
  FLA_Error e_val;
  dim_t     m_a;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );  

  e_val = FLA_Check_floating_object( a );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( a );
  FLA_Check_error_code( e_val );

  m_a = FLA_Obj_vector_dim( a );
  
  if ( m_a > 1 ) 
  {
    e_val = FLA_Check_floating_object( b );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_if_vector( b );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_identical_object_datatype( a, d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( a, e );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( e );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( d, m_a );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( e, m_a );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

