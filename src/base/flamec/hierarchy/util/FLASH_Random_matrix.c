
#include "FLAME.h"

FLA_Error FLASH_Random_matrix( FLA_Obj H )
{
	FLA_Obj F;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Randomize the flat matrix object.
	FLA_Random_matrix( F );

	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( F, H );

	// Free the temporary flat object.
	FLA_Obj_free( &F );

	return FLA_SUCCESS;
}

