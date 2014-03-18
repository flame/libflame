
#include "FLAME.h"

fla_scalr_t* fla_scalr_cntl_blas;

void FLA_Scalr_cntl_init()
{
	// Create a control tree that assumes A is small.
	fla_scalr_cntl_blas  = FLA_Cntl_scalr_obj_create( FLA_FLAT,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );
}

void FLA_Scalr_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_scalr_cntl_blas );
}

