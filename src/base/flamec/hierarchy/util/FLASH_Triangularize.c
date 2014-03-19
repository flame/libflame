/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
	FLA_Error r_val;
	FLA_Obj   A_flat;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );

	// Triangularize the flat matrix object as specified by uplo and diag.
	r_val = FLA_Triangularize( uplo, diag, A_flat );
	
	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( A_flat, A );

	// Free the temporary flat object.
	FLA_Obj_free( &A_flat );
	
	return r_val;
}

