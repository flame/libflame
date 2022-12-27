/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_UDdate_UT_inc_create_hier_matrices( FLA_Obj R_flat, FLA_Obj C_flat, FLA_Obj D_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* R, FLA_Obj* C, FLA_Obj* D, FLA_Obj* T, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        m_T, n_T;
	dim_t        m_W, n_W;
	dim_t        m_C;
	dim_t        m_D;
	
	// *** The current UDdate_UT_inc algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1. We check for that here
	// because we anticipate that we'll use a more general algorithm in the
	// future, and we don't want to forget to remove the constraint. ***
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_UDdate_UT_inc() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Create hierarchical copy of matrices R_flat, C_flat, and D_flat.
	FLASH_Obj_create_hier_copy_of_flat( R_flat, depth, b_flash, R );
	FLASH_Obj_create_hier_copy_of_flat( C_flat, depth, b_flash, C );
	FLASH_Obj_create_hier_copy_of_flat( D_flat, depth, b_flash, D );

	// Query the datatype of matrix R_flat.
	datatype = FLA_Obj_datatype( R_flat );
	
	// If the user passed in zero for b_alg, then we need to set the
	// algorithmic (inner) blocksize to a reasonable default value.
	if ( b_alg == 0 )
	{
		b_alg = FLASH_UDdate_UT_inc_determine_alg_blocksize( *R );
	}

	// Determine the element (not scalar) dimensions of the new hierarchical
	// matrix T. By using the element dimensions, we will probably allocate
	// more storage than we actually need (at the bottom and right edge cases)
	// but this is simpler than computing the exact amount and the excess
	// storage is usually small in practice.
	n_T = FLA_Obj_width( *R );
	m_C = FLA_Obj_length( *C );
	m_D = FLA_Obj_length( *D );
	m_T = fla_max( m_C, m_D );

	// Create hierarchical matrix T, with element dimensions conformal to the
	// the larger of C and D, where each block is b_alg-by-b_flash.
	FLASH_Obj_create_ext( datatype, m_T * b_alg, n_T * b_flash[0], 
	                      depth, &b_alg, b_flash, 
	                      T );

	// Determine the element (not scalar) dimensions of the new hierarchical
	// matrix W. The element length and width will be identical to that of R.
	// Once again, we will probably allocate excess storage, but we consider
	// this to be small.
	m_W = FLA_Obj_length( *R );
	n_W = FLA_Obj_width( *R );
	   
	// Create hierarchical matrix W, with element dimensions conformal to R,
	// where each block is b_alg-by-b_flash.
	FLASH_Obj_create_ext( datatype, m_W * b_alg, n_W * b_flash[0], 
	                      depth, &b_alg, b_flash, 
	                      W );

	return FLA_SUCCESS;
}


dim_t FLASH_UDdate_UT_inc_determine_alg_blocksize( FLA_Obj R )
{
	dim_t b_alg;
	dim_t b_flash;

	// Acquire the storage blocksize.
	b_flash = FLA_Obj_length( *FLASH_OBJ_PTR_AT( R ) );

	// Scale the storage blocksize by a pre-defined scalar to arrive at a
	// reasonable algorithmic blocksize, but make sure it's at least 1.
	b_alg = ( dim_t ) fla_max( ( double ) b_flash * FLA_UDDATE_INNER_TO_OUTER_B_RATIO, 1 );

	return b_alg;
}

