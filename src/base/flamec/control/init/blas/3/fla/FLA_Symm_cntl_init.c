/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* fla_scal_cntl_blas;
extern fla_gemm_t* fla_gemm_cntl_blas;

fla_symm_t*        fla_symm_cntl_blas = NULL;
fla_symm_t*        fla_symm_cntl_bp = NULL;
fla_symm_t*        fla_symm_cntl_mp = NULL;
fla_symm_t*        fla_symm_cntl_mm = NULL;
fla_blocksize_t*   fla_symm_var1_bsize = NULL;
fla_blocksize_t*   fla_symm_var9_bsize = NULL;

void FLA_Symm_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_symm_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_symm_var9_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A and B are b x b blocks.
	fla_symm_cntl_blas  = FLA_Cntl_symm_obj_create( FLA_FLAT, 
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	fla_symm_cntl_bp    = FLA_Cntl_symm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT9,
	                                                fla_symm_var9_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_symm_cntl_blas,
	                                                NULL,
	                                                NULL );

	// Create a control tree that assumes A is large and B is a panel.
	fla_symm_cntl_mp    = FLA_Cntl_symm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT1,
	                                                fla_symm_var1_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_symm_cntl_blas,
	                                                fla_gemm_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create a control tree that assumes A and B are both large.
	fla_symm_cntl_mm    = FLA_Cntl_symm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT9,
	                                                fla_symm_var9_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_symm_cntl_mp,
	                                                NULL,
	                                                NULL );

}

void FLA_Symm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_symm_cntl_blas );

	FLA_Cntl_obj_free( fla_symm_cntl_bp );
	FLA_Cntl_obj_free( fla_symm_cntl_mp );
	FLA_Cntl_obj_free( fla_symm_cntl_mm );

	FLA_Blocksize_free( fla_symm_var1_bsize );
	FLA_Blocksize_free( fla_symm_var9_bsize );
}

