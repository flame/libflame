
#include "FLAME.h"


FLA_Error FLA_Bsvd_ext( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H,
                        FLA_Svd_type jobu, FLA_Obj U,
                        FLA_Svd_type jobv, FLA_Obj V,
                        FLA_Bool apply_Uh2C, FLA_Obj C )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     b_alg      = 512;
  dim_t     m_d        = FLA_Obj_vector_dim( d );
  FLA_Obj   W;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bsvd_ext_check( uplo, d, e, G, H, jobu, U, jobv, V, apply_Uh2C, C );

  // Partition U and V properly to account for the diagonal dimension.
  if ( jobu == FLA_SVD_VECTORS_MIN_COPY || 
       jobu == FLA_SVD_VECTORS_MIN_OVERWRITE )
    FLA_Part_1x2( U, &U, &W, m_d, FLA_LEFT );
  
  if ( jobv == FLA_SVD_VECTORS_MIN_COPY || 
       jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
    FLA_Part_1x2( V, &V, &W, m_d, FLA_LEFT );

  // when uplo is lower triangular, swap U and V.
  if ( uplo == FLA_LOWER_TRIANGULAR )
    exchange( U, V, W );
  
  // [ U1 d e V1 ] = Bidiag( A );
  // [ U2 d V2 ] = Bsvd( diag( d ) + subdiag( uplo, e ) ); 
  // Then, 
  // A = U1 U2 d V2^H V1^H;
  // solving AX = C, it would be convenient that
  // d (V1 V2)^H = apply(U2^H) apply(U1^H) C
  // Here apply(.) is implicitly made and (V1 V2)^H is explicitly overwritten on the A.
  // So, apply_Uh2C is a memory efficient way to solve a system of equations.
  r_val = FLA_Bsvd_ext_opt_var1( n_iter_max, 
                                 d, e, G, H, 
                                 jobu, U, jobv, V, 
                                 apply_Uh2C, C,
                                 b_alg );

  return r_val;
}

