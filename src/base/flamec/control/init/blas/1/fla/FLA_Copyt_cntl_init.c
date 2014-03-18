
#include "FLAME.h"

fla_copyt_t* fla_copyt_cntl_blas;

void FLA_Copyt_cntl_init()
{
	// Create a control tree that assumes A and B are small.
	fla_copyt_cntl_blas = FLA_Cntl_copyt_obj_create( FLA_FLAT,
	                                                 FLA_SUBPROBLEM,
	                                                 NULL,
	                                                 NULL );
}

void FLA_Copyt_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_copyt_cntl_blas );
}

