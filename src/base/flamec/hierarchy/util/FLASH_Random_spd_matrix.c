/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj H )
{
	FLA_Obj   F;
	FLA_Error e_val;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_valid_uplo( uplo );
		FLA_Check_error_code( e_val );
	}

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Randomize the flat matrix object to be SPD.
	FLA_Random_spd_matrix( uplo, F );

	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( F, H );

	// Free the temporary flat object.
	FLASH_Obj_free( &F );

	return FLA_SUCCESS;
}

