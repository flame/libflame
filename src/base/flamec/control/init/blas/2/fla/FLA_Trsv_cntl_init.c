
#include "FLAME.h"

fla_trsv_t* fla_trsv_cntl_blas;

void FLA_Trsv_cntl_init()
{
	// Create a control tree that assumes A is small.
	fla_trsv_cntl_blas  = FLA_Cntl_trsv_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL );
}

void FLA_Trsv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_trsv_cntl_blas );
}

