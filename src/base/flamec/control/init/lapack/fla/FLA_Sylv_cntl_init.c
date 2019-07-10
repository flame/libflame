/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t* fla_gemm_cntl_blas;

fla_sylv_t*        fla_sylv_cntl_leaf = NULL;
fla_sylv_t*        fla_sylv_cntl_mb = NULL;
fla_sylv_t*        fla_sylv_cntl = NULL;
fla_blocksize_t*   fla_sylv_bsize = NULL;

void FLA_Sylv_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_sylv_bsize       = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke LAPACK.
	fla_sylv_cntl_leaf   = FLA_Cntl_sylv_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                 FLA_BLOCKED_EXTERN, 
#else
	                                                 FLA_UNB_OPT_VARIANT1,
#endif
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 15.
	fla_sylv_cntl_mb     = FLA_Cntl_sylv_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT15,
	                                                 fla_sylv_bsize,
	                                                 fla_sylv_cntl_leaf,
	                                                 NULL,
	                                                 NULL,
	                                                 fla_gemm_cntl_blas,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 17.
	fla_sylv_cntl        = FLA_Cntl_sylv_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT17,
	                                                 fla_sylv_bsize,
	                                                 fla_sylv_cntl_mb,
	                                                 NULL,
	                                                 NULL,
	                                                 fla_gemm_cntl_blas,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );
}

void FLA_Sylv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_sylv_cntl_leaf );
	FLA_Cntl_obj_free( fla_sylv_cntl_mb );
	FLA_Cntl_obj_free( fla_sylv_cntl );

	FLA_Blocksize_free( fla_sylv_bsize );
}

