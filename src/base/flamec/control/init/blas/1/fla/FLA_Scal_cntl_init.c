
#include "FLAME.h"

fla_scal_t* fla_scal_cntl_blas;

void FLA_Scal_cntl_init()
{
	// Create a control tree that assumes A is small.
	fla_scal_cntl_blas  = FLA_Cntl_scal_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL );
}

void FLA_Scal_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_scal_cntl_blas );
}

