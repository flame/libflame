/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t* flash_gemm_cntl_op_bp;
extern fla_trsm_t* flash_trsm_cntl_bp;

fla_trinv_t*       flash_trinv_cntl_leaf = NULL;
fla_trinv_t*       flash_trinv_cntl = NULL;
fla_blocksize_t*   flash_trinv_bsize = NULL;

void FLASH_Trinv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_trinv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_trinv_cntl_leaf   = FLA_Cntl_trinv_obj_create( FLA_HIER,
	                                                     FLA_SUBPROBLEM,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL );

	// Create a control tree that assumes A is large.
	flash_trinv_cntl        = FLA_Cntl_trinv_obj_create( FLA_HIER,
	                                                     FLA_BLOCKED_VARIANT3, 
	                                                     flash_trinv_bsize,
	                                                     flash_trinv_cntl_leaf,
	                                                     NULL,
	                                                     flash_trsm_cntl_bp,
	                                                     flash_trsm_cntl_bp,
	                                                     flash_gemm_cntl_op_bp );
}

void FLASH_Trinv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_trinv_cntl_leaf );
	FLA_Cntl_obj_free( flash_trinv_cntl );

	FLA_Blocksize_free( flash_trinv_bsize );
}

