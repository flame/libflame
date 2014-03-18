
#include "FLAME.h"

FLA_Error FLASH_Norm1( FLA_Obj H, FLA_Obj norm )
{
	FLA_Obj F;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) )
	{
		FLA_Set( FLA_ZERO, norm );
		return FLA_SUCCESS;
	}

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Compute the 1-norm of F and store it in norm.
	FLA_Norm1( F, norm );
	
	// Free the temporary flat object.
	FLA_Obj_free( &F );
	
	return FLA_SUCCESS;
}

