/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

// jobu,   jobv   [in] indicates operations for U and V.
// transu, transv [in] indicates shape of output matrices (e.g.,  U or V^H )
// Dimensions of U and V should match to corresponding trans parameters.
//
// For job == FLA_SVD_VECTORS_MIN_OVERWRITE, trans must be set properly:
//
//   jobu = FLA_SVD_VECTORS_MIN_OVERWRITE, transu = FLA_NO_TRANSPOSE;
//   jobv = FLA_SVD_VECTORS_MIN_OVERWRITE, transv = FLA_CONJ_TRANSPOSE;
//
// Also note that jobu(v) == FLA_SVD_VECTORS_NONE does not reference U (or V).
// This implies that it is illegal to access their base objects and associated 
// data (cs, rs, and datatype).
FLA_Error FLA_Svd_ext( FLA_Svd_type jobu, FLA_Trans transu, 
                       FLA_Svd_type jobv, FLA_Trans transv, 
                       FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     k_accum    = 32;
  dim_t     b_alg      = 512;
  dim_t     min_m_n    = FLA_Obj_min_dim( A );
  dim_t     m_A        = FLA_Obj_length( A );
  dim_t     n_A        = FLA_Obj_width( A );
  FLA_Bool  u_flipped  = FALSE;
  FLA_Bool  v_flipped  = FALSE;
  FLA_Obj   W;           // Dummy variable for partitioning of matrices.

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Svd_ext_check( jobu, transu, jobv, transv, A, s, U, V );

  // Transpose U and V to match dimensions used in SVD. 
  if ( ( transu == FLA_TRANSPOSE        || transu == FLA_CONJ_TRANSPOSE            ) &&
       ( jobu   != FLA_SVD_VECTORS_NONE && jobu   != FLA_SVD_VECTORS_MIN_OVERWRITE ) )
  {
    FLA_Obj_flip_base( &U );
    FLA_Obj_flip_view( &U );
    u_flipped = TRUE;
  }
  if ( ( transv == FLA_TRANSPOSE        || transv == FLA_CONJ_TRANSPOSE            ) &&
       ( jobv   != FLA_SVD_VECTORS_NONE && jobv   != FLA_SVD_VECTORS_MIN_OVERWRITE ) )
  {
    FLA_Obj_flip_base( &V );
    FLA_Obj_flip_view( &V );
    v_flipped = TRUE;
  }

  // Partition U and V if necessary.
  if ( jobu == FLA_SVD_VECTORS_MIN_COPY ) FLA_Part_1x2( U, &U, &W, min_m_n, FLA_LEFT );
  if ( jobv == FLA_SVD_VECTORS_MIN_COPY ) FLA_Part_1x2( V, &V, &W, min_m_n, FLA_LEFT );

  if ( m_A >= n_A )
  {
    r_val = FLA_Svd_ext_u_unb_var1( jobu, jobv, 
                                    n_iter_max, 
                                    A, s, U, V, 
                                    k_accum, b_alg );

    // Recover U and V.
    if ( u_flipped == TRUE ) 
    {
      if ( FLA_Obj_is_complex( U ) )
        FLA_Conjugate( U );
      FLA_Obj_flip_base( &U );
    }
    if ( v_flipped == TRUE ) 
    {
      if ( FLA_Obj_is_complex( V ) )
        FLA_Conjugate( V );
      FLA_Obj_flip_base( &V );
    }
  }
  else
  {
    // Flip A and exchange U and V parameters.
    FLA_Obj_flip_base( &A );
    FLA_Obj_flip_view( &A );

    // Note that U and V are also swapped.
    r_val = FLA_Svd_ext_u_unb_var1( jobv, jobu, 
                                    n_iter_max, 
                                    A, s, V, U, 
                                    k_accum, b_alg );
    
    // Recover A.
    FLA_Obj_flip_base( &A );

    // Recover U and V. Consider a case that U and V are not created. 
    if ( u_flipped == TRUE ) 
      FLA_Obj_flip_base( &U );
    else if ( jobu != FLA_SVD_VECTORS_NONE && 
              jobu != FLA_SVD_VECTORS_MIN_OVERWRITE ) 
      if ( FLA_Obj_is_complex( U ) )
        FLA_Conjugate( U );
    
    if ( v_flipped == TRUE ) 
      FLA_Obj_flip_base( &V );
    else if ( jobv != FLA_SVD_VECTORS_NONE && 
              jobv != FLA_SVD_VECTORS_MIN_OVERWRITE )
      if ( FLA_Obj_is_complex( V ) )
        FLA_Conjugate( V );
  }

  return r_val;
}

