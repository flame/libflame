
#include "FLAME.h"

FLA_Error FLASH_Shift_diag( FLA_Conj conj, FLA_Obj sigma, FLA_Obj H )
{
	FLA_Obj F;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Shift the diagonal of the flat matrix object by sigma.
	FLA_Shift_diag( conj, sigma, F );

	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( F, H );
	
	// Free the temporary flat object.
	FLA_Obj_free( &F );
	
	return FLA_SUCCESS;
}

