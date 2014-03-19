/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_check( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e,
                          FLA_Obj G, FLA_Obj H,
                          FLA_Svd_type jobu, FLA_Obj U,
                          FLA_Svd_type jobv, FLA_Obj V )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     m_d, m_e;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( d );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( d );
  FLA_Check_error_code( e_val );

  m_d = FLA_Obj_vector_dim( d );
  m_e = ( m_d - 1 );

  if ( m_e > 0 )
  {
    e_val = FLA_Check_real_object( e );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_nonconstant_object( e );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_identical_object_datatype( d, e );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_vector_dim( e, m_e );
    FLA_Check_error_code( e_val );
  }
  
  if ( m_e > 0 )
  {
    e_val = FLA_Check_complex_object( G );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_identical_object_precision( d, G );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_object_length_equals( G, m_e );    
    FLA_Check_error_code( e_val );
  }

  if ( m_e > 0 )
  {
    e_val = FLA_Check_complex_object( H );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_identical_object_precision( d, H );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_length_equals( H, m_e );    
    FLA_Check_error_code( e_val );
  }

  FLA_Check_valid_svd_type( jobu );
  FLA_Check_error_code( e_val );

  FLA_Check_valid_svd_type( jobv );
  FLA_Check_error_code( e_val );

  if ( jobu != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( U ) == FALSE )
  {
    FLA_Check_identical_object_precision( d, U );
    FLA_Check_error_code( e_val );
  }
  
  if ( jobv != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( V ) == FALSE )
  {
    FLA_Check_identical_object_precision( d, V );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

