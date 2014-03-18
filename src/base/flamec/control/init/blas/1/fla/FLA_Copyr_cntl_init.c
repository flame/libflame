
#include "FLAME.h"

fla_copyr_t* fla_copyr_cntl_blas;

void FLA_Copyr_cntl_init()
{
	// Create a control tree that assumes A and B are small.
	fla_copyr_cntl_blas = FLA_Cntl_copyr_obj_create( FLA_FLAT,
	                                                 FLA_SUBPROBLEM,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );
}

void FLA_Copyr_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_copyr_cntl_blas );
}

