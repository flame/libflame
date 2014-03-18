
#include "FLAME.h"

FLA_Error FLASH_Hermitianize( FLA_Uplo uplo, FLA_Obj A )
{
	FLA_Error r_val;
	FLA_Obj   A_flat;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );

	// Hermitianize the flat matrix object as specified by uplo.
	r_val = FLA_Hermitianize( uplo, A_flat );
	
	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( A_flat, A );

	// Free the temporary flat object.
	FLA_Obj_free( &A_flat );
	
	return r_val;
}

