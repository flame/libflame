/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_herk_t* flash_herk_cntl_op;
extern fla_trsm_t* flash_trsm_cntl_bp;

fla_chol_t*        flash_chol_cntl_leaf = NULL;
fla_chol_t*        flash_chol_cntl = NULL;
fla_blocksize_t*   flash_chol_bsize = NULL;

void FLASH_Chol_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_chol_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_chol_cntl_leaf   = FLA_Cntl_chol_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is large.
	flash_chol_cntl        = FLA_Cntl_chol_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3, 
	                                                   flash_chol_bsize,
	                                                   flash_chol_cntl_leaf,
	                                                   flash_herk_cntl_op,
	                                                   flash_trsm_cntl_bp,
	                                                   NULL );
}

void FLASH_Chol_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_chol_cntl_leaf );
	FLA_Cntl_obj_free( flash_chol_cntl );

	FLA_Blocksize_free( flash_chol_bsize );
}

