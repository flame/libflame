
#include "FLAME.h"

fla_copy_t* fla_copy_cntl_blas;

void FLA_Copy_cntl_init()
{
	// Create a control tree that assumes A and B are small.
	fla_copy_cntl_blas  = FLA_Cntl_copy_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL );
}

void FLA_Copy_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_copy_cntl_blas );
}

