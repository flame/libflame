
#include "FLAME.h"

fla_axpyt_t* fla_axpyt_cntl_blas;

void FLA_Axpyt_cntl_init()
{
	// Create a control tree that assumes A and B are small.
	fla_axpyt_cntl_blas  = FLA_Cntl_axpyt_obj_create( FLA_FLAT,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL );
}

void FLA_Axpyt_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_axpyt_cntl_blas );
}

