
#include "FLAME.h"

fla_gemv_t* fla_gemv_cntl_blas;

void FLA_Gemv_cntl_init()
{
	// Create a control tree that assumes A is small.
	fla_gemv_cntl_blas  = FLA_Cntl_gemv_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL );
}

void FLA_Gemv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_gemv_cntl_blas );
}

