
#include "FLAME.h"

fla_axpy_t* fla_axpy_cntl_blas;

void FLA_Axpy_cntl_init()
{
	// Create a control tree that assumes A and B are small.
	fla_axpy_cntl_blas  = FLA_Cntl_axpy_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL );
}

void FLA_Axpy_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_axpy_cntl_blas );
}

