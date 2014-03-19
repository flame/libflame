/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_LQ_UT_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, FLA_Obj* A, FLA_Obj* TW )
{
	FLA_Datatype datatype;
	dim_t        m, n;
	dim_t        min_m_n;
	
	// *** The current LQ_UT algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1. We check for that here
	// because we anticipate that we'll use a more general algorithm in the
	// future, and we don't want to forget to remove the constraint. ***
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_LQ_UT() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Create hierarchical copy of matrix A_flat.
	FLASH_Obj_create_hier_copy_of_flat( A_flat, depth, b_flash, A );

	// Query the datatype of matrix A_flat.
	datatype = FLA_Obj_datatype( A_flat );
	
	// Query the minimum dimension of A_flat.
	min_m_n = FLA_Obj_min_dim( A_flat );

	// Set the m and n dimensions of TW to be min_m_n.
	m = min_m_n;
	n = min_m_n;

	// Create hierarchical matrices T and W.
	FLASH_Obj_create_ext( datatype, m, n, 
	                      depth, b_flash, b_flash, 
	                      TW );
	   
	return FLA_SUCCESS;
}

