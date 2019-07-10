/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_scalr_t* fla_scalr_cntl_blas = NULL;

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

