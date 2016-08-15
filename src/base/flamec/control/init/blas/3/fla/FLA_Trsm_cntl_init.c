/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* fla_scal_cntl_blas;
extern fla_gemm_t* fla_gemm_cntl_blas;

fla_trsm_t*        fla_trsm_cntl_blas = NULL;
fla_trsm_t*        fla_trsm_cntl_bp = NULL;
fla_trsm_t*        fla_trsm_cntl_mp = NULL;
fla_trsm_t*        fla_trsm_cntl_mm = NULL;
fla_blocksize_t*   fla_trsm_var2_bsize = NULL;
fla_blocksize_t*   fla_trsm_var3_bsize = NULL;

void FLA_Trsm_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_trsm_var2_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_trsm_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A and B are b x b blocks.
	fla_trsm_cntl_blas  = FLA_Cntl_trsm_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	fla_trsm_cntl_bp    = FLA_Cntl_trsm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_trsm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trsm_cntl_blas,
	                                                NULL );

	// Create a control tree that assumes A is large and B is a panel.
	fla_trsm_cntl_mp    = FLA_Cntl_trsm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT2,
	                                                fla_trsm_var2_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trsm_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create a control tree that assumes A and B are both large.
	fla_trsm_cntl_mm    = FLA_Cntl_trsm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_trsm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trsm_cntl_mp,
	                                                NULL );
}

void FLA_Trsm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_trsm_cntl_blas );

	FLA_Cntl_obj_free( fla_trsm_cntl_bp );
	FLA_Cntl_obj_free( fla_trsm_cntl_mp );
	FLA_Cntl_obj_free( fla_trsm_cntl_mm );

	FLA_Blocksize_free( fla_trsm_var2_bsize );
	FLA_Blocksize_free( fla_trsm_var3_bsize );
}

