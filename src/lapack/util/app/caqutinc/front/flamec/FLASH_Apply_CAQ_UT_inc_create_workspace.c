/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Apply_CAQ_UT_inc_create_workspace( dim_t p, FLA_Obj TW, FLA_Obj B, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        depth;
	dim_t        b_alg;
	dim_t        b_flash;
	dim_t        m, n;

	// Query the depth.
	depth = FLASH_Obj_depth( TW );
	
	// *** The current Apply_CAQ_UT_inc algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1.
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_Apply_CAQ_UT_inc() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Query the datatype of matrix TW.
	datatype = FLA_Obj_datatype( TW );
	
	// Inspect the length of a the top-left element of TW to get the
	// algorithmic blocksize we'll use throughout the Apply_CAQ_UT_inc
	// algorithm.
	b_alg = FLASH_Obj_scalar_length_tl( TW );

	// The width of the top-left element gives us the storage blocksize.
	b_flash = FLASH_Obj_scalar_width_tl( TW );

	// The element length of W need to be p: one panel for each
    // factorized subproblem.
	m = p;

	// Query the element (not scalar) width of the right-hand side
	// matrix B. This is done so we can create W with full blocks for the
	// right "edge cases" of B.
	n = FLA_Obj_width( B );

	// Create hierarchical matrix W.
	FLASH_Obj_create_ext( datatype, m * b_alg, n * b_flash, 
	                      depth, &b_alg, &b_flash, 
	                      W );
	   
	return FLA_SUCCESS;
}

