/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t* fla_gemm_cntl_blas;
extern fla_herk_t* fla_herk_cntl_blas;
extern fla_trmm_t* fla_trmm_cntl_blas;

fla_ttmm_t*        fla_ttmm_cntl_leaf = NULL;
fla_ttmm_t*        fla_ttmm_cntl = NULL;
fla_blocksize_t*   fla_ttmm_var1_bsize = NULL;

void FLA_Ttmm_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_ttmm_var1_bsize  = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke LAPACK.
	fla_ttmm_cntl_leaf   = FLA_Cntl_ttmm_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                 FLA_BLOCKED_EXTERN, 
#else
	                                                 FLA_UNB_OPT_VARIANT2,
#endif
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 1.
	fla_ttmm_cntl        = FLA_Cntl_ttmm_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT1, 
	                                                 fla_ttmm_var1_bsize,
	                                                 fla_ttmm_cntl_leaf,
	                                                 fla_herk_cntl_blas,
	                                                 fla_trmm_cntl_blas,
	                                                 fla_gemm_cntl_blas );
}

void FLA_Ttmm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_ttmm_cntl_leaf );
	FLA_Cntl_obj_free( fla_ttmm_cntl );

	FLA_Blocksize_free( fla_ttmm_var1_bsize );
}

