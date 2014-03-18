
#include "FLAME.h"

FLA_Error FLA_Svd( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     k_accum    = 32;
  dim_t     b_alg      = 512;
  dim_t     min_m_n    = FLA_Obj_min_dim( A );
  dim_t     m_A        = FLA_Obj_length( A );
  dim_t     n_A        = FLA_Obj_width( A );
  FLA_Obj   W;         // Dummy variable for partitioning of matrices.

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Svd_check( jobu, jobv, A, s, U, V );

  // Partition U and V if necessary.
  if ( jobu == FLA_SVD_VECTORS_MIN_COPY ) FLA_Part_1x2( U, &U, &W, min_m_n, FLA_LEFT );
  if ( jobv == FLA_SVD_VECTORS_MIN_COPY ) FLA_Part_1x2( V, &V, &W, min_m_n, FLA_LEFT );

  // Use extension version
  if ( m_A >= n_A )
  {
    r_val = FLA_Svd_ext_u_unb_var1( jobu, jobv, 
                                    n_iter_max, 
                                    A, s, U, V,
                                    k_accum, b_alg );
  }
  else
  {
    // Flip A and change U and V
    FLA_Obj_flip_base( &A );
    FLA_Obj_flip_view( &A );
    
    r_val = FLA_Svd_ext_u_unb_var1( jobu, jobv, 
                                    n_iter_max, 
                                    A, s, V, U,
                                    k_accum, b_alg );
    
    // Recover A and conjugate U and V for complex cases
    FLA_Obj_flip_base( &A );
    
    if ( FLA_Obj_is_complex( A ) )
    {
      if ( jobu != FLA_SVD_VECTORS_NONE ) FLA_Conjugate( U );
      if ( jobv != FLA_SVD_VECTORS_NONE ) FLA_Conjugate( V );
    }
  }

  return r_val;
}
