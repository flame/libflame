/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, 
                    FLA_Obj G, FLA_Obj H, 
                    FLA_Svd_type jobu, FLA_Obj U, 
                    FLA_Svd_type jobv, FLA_Obj V )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     b_alg      = 512;
  dim_t     m_d        = FLA_Obj_vector_dim( d );
  FLA_Obj   C; // dummy variable

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bsvd_check( uplo, d, e, G, H, jobu, U, jobv, V );

  // Partition U and V properly to account for the diagonal dimension.
  if ( jobu == FLA_SVD_VECTORS_MIN_COPY || 
       jobu == FLA_SVD_VECTORS_MIN_OVERWRITE )
    FLA_Part_1x2( U, &U, &C, m_d, FLA_LEFT );
  
  if ( jobv == FLA_SVD_VECTORS_MIN_COPY || 
       jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
    FLA_Part_1x2( V, &V, &C, m_d, FLA_LEFT );

  // The current Bsvd is implemented for the upper triangular matrix only.
  // If uplo is lower triangular, swap U and V.
  if ( uplo == FLA_LOWER_TRIANGULAR )
    exchange( U, V, C );
  
  r_val = FLA_Bsvd_ext_opt_var1 ( n_iter_max, d, e, G, H, 
                                  jobu, U, jobv, V, 
                                  FALSE, C,
                                  b_alg );

  return r_val;
}

