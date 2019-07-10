/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t* flash_gemm_cntl_pm_bp;
extern fla_gemm_t* flash_gemm_cntl_ip_bb;

fla_sylv_t*        flash_sylv_cntl_leaf = NULL;
fla_sylv_t*        flash_sylv_cntl_mb = NULL;
fla_sylv_t*        flash_sylv_cntl = NULL;
fla_blocksize_t*   flash_sylv_bsize = NULL;

void FLASH_Sylv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_sylv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_sylv_cntl_leaf   = FLA_Cntl_sylv_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
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

	// Create a control tree that assumes A is a matrix and B is a block.
	flash_sylv_cntl_mb     = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT17,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_leaf,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_ip_bb,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is a matrix and B is a matrix.
	flash_sylv_cntl        = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT15,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_mb,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_pm_bp,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );
}

void FLASH_Sylv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_sylv_cntl_leaf );
	FLA_Cntl_obj_free( flash_sylv_cntl_mb );
	FLA_Cntl_obj_free( flash_sylv_cntl );

	FLA_Blocksize_free( flash_sylv_bsize );
}

