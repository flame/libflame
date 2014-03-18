
#include "FLAME.h"

FLA_Error FLA_Svd_ext_check( FLA_Svd_type jobu, FLA_Trans transu, 
                             FLA_Svd_type jobv, FLA_Trans transv, 
                             FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_svd_type( jobu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_svd_type( jobv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_svd_type_combination( jobu, jobv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( transu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( transv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_svd_type_and_trans_combination( jobu, transu, jobv, transv );
  FLA_Check_error_code( e_val );

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

  if ( jobu != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( U ) == FALSE ) 
  {
    e_val = FLA_Check_identical_object_datatype( A, U );
    FLA_Check_error_code( e_val );

    if ( transu == FLA_NO_TRANSPOSE || transu == FLA_CONJ_NO_TRANSPOSE )
    {
      e_val = FLA_Check_object_length_equals( U, FLA_Obj_length( A ) );
      FLA_Check_error_code( e_val );
    }
    else
    {
      e_val = FLA_Check_object_width_equals( U, FLA_Obj_length( A ) );
      FLA_Check_error_code( e_val );
    }
  }
  
  if ( jobv != FLA_SVD_VECTORS_NONE && FLA_Obj_has_zero_dim( V ) == FALSE ) 
  {
    e_val = FLA_Check_identical_object_datatype( A, V );
    FLA_Check_error_code( e_val );

    if ( transv == FLA_NO_TRANSPOSE || transv == FLA_CONJ_NO_TRANSPOSE )
    {
      e_val = FLA_Check_object_length_equals( V, FLA_Obj_width( A ) );
      FLA_Check_error_code( e_val );
    }
    else
    {
      e_val = FLA_Check_object_width_equals( V, FLA_Obj_width( A ) );
      FLA_Check_error_code( e_val );
    }
  }

  return FLA_SUCCESS;
}

