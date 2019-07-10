/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t*  flash_gemm_cntl_bp_bb;
extern fla_trsm_t*  flash_trsm_cntl_bp;
extern fla_appiv_t* flash_appiv_cntl_bp;

fla_lu_t*           flash_lu_incpiv_cntl_leaf = NULL;
fla_lu_t*           flash_lu_incpiv_cntl = NULL;
fla_blocksize_t*    flash_lu_incpiv_bsize = NULL;

void FLASH_LU_incpiv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_lu_incpiv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_lu_incpiv_cntl_leaf   = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                      FLA_SUBPROBLEM,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL );

	// Create a control tree that assumes A is large.
	flash_lu_incpiv_cntl        = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      flash_lu_incpiv_bsize,
	                                                      flash_lu_incpiv_cntl_leaf,
	                                                      flash_gemm_cntl_bp_bb,
	                                                      NULL,
	                                                      NULL,
	                                                      flash_trsm_cntl_bp,
	                                                      flash_trsm_cntl_bp,
	                                                      flash_appiv_cntl_bp,
	                                                      NULL );
}

void FLASH_LU_incpiv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_lu_incpiv_cntl_leaf );
	FLA_Cntl_obj_free( flash_lu_incpiv_cntl );

	FLA_Blocksize_free( flash_lu_incpiv_bsize );
}

