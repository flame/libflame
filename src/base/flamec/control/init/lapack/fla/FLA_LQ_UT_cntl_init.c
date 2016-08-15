/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_apqut_t* fla_apqut_cntl_leaf;

fla_lqut_t*         fla_lqut_cntl_unb = NULL;
fla_lqut_t*         fla_lqut_cntl_leaf = NULL;

fla_blocksize_t*    fla_lqut_var1_bsize_leaf = NULL;

void FLA_LQ_UT_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_lqut_var1_bsize_leaf = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_lqut_var1_bsize_leaf, FLA_LQ_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke unblocked variant 2.
	fla_lqut_cntl_unb  = FLA_Cntl_lqut_obj_create( FLA_FLAT,
	                                               FLA_UNB_OPT_VARIANT2,
	                                               NULL,
	                                               NULL,
	                                               NULL );

	// Create a control tree for small-to-medium sequential problems and
	// as the means to compute on FLASH blocks.
	fla_lqut_cntl_leaf = FLA_Cntl_lqut_obj_create( FLA_FLAT, 
	                                               FLA_BLOCKED_VARIANT1,
	                                               fla_lqut_var1_bsize_leaf,
	                                               fla_lqut_cntl_unb,
	                                               fla_apqut_cntl_leaf );
}

void FLA_LQ_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_lqut_cntl_unb );
	FLA_Cntl_obj_free( fla_lqut_cntl_leaf );

	FLA_Blocksize_free( fla_lqut_var1_bsize_leaf );
}

