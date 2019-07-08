/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLA_Gemv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	// Create a control tree that assumes A is small.
	FLA_cntl_flamec_init_i->fla_gemv_cntl_blas  = FLA_Cntl_gemv_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL );
}

void FLA_Gemv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flamec_init_i->fla_gemv_cntl_blas );
}
#endif

fla_gemv_t* fla_gemv_cntl_blas = NULL;

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

