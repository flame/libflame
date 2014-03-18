
#include "FLAME.h"

FLA_Error FLASH_Set( FLA_Obj alpha, FLA_Obj H )
{
	FLA_Obj F;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Scale the flat matrix object by alpha.
	FLA_Set( alpha, F );

	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( F, H );
	
	// Free the temporary flat object.
	FLA_Obj_free( &F );
	
	return FLA_SUCCESS;
}

