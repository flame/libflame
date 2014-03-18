
#include "FLAME.h"

double FLASH_Max_elemwise_diff( FLA_Obj A, FLA_Obj B )
{
	FLA_Obj A_flat, B_flat;
	double  max_diff;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( A ) ) return -1.0;

	// Create a temporary flat copy of the hierarchical objects.
	FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
	FLASH_Obj_create_flat_copy_of_hier( B, &B_flat );

	// Get the maximum element-wise diff.
	max_diff = FLA_Max_elemwise_diff( A_flat, B_flat );
	
	// Free the temporary flat objects.
	FLA_Obj_free( &A_flat );
	FLA_Obj_free( &B_flat );
	
	return max_diff;
}

