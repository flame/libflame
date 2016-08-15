/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_apqut_t*  flash_apqut_cntl_blas;

fla_lqut_t*          flash_lqut_cntl_leaf = NULL;
fla_lqut_t*          flash_lqut_cntl = NULL;

fla_blocksize_t*     flash_lqut_var3_bsize = NULL;

void FLASH_LQ_UT_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_lqut_var3_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to compute the subproblem.
	flash_lqut_cntl_leaf = FLA_Cntl_lqut_obj_create( FLA_HIER,
	                                                 FLA_SUBPROBLEM, 
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 2.
	flash_lqut_cntl      = FLA_Cntl_lqut_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT3,
	                                                 flash_lqut_var3_bsize,
	                                                 flash_lqut_cntl_leaf,
	                                                 flash_apqut_cntl_blas );
}

void FLASH_LQ_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_lqut_cntl_leaf );
	FLA_Cntl_obj_free( flash_lqut_cntl );

	FLA_Blocksize_free( flash_lqut_var3_bsize );
}

