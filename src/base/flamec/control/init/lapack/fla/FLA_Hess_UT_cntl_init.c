/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_hessut_t*       fla_hessut_cntl_leaf = NULL;

fla_blocksize_t*    fla_hessut_bsize_leaf = NULL;

void FLA_Hess_UT_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_hessut_bsize_leaf = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_hessut_bsize_leaf, FLA_HESS_INNER_TO_OUTER_B_RATIO );

	// Create a control tree for small-to-medium sequential problems.
	fla_hessut_cntl_leaf = FLA_Cntl_hessut_obj_create( FLA_FLAT, 
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_hessut_bsize_leaf );
}

void FLA_Hess_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_hessut_cntl_leaf );

	FLA_Blocksize_free( fla_hessut_bsize_leaf );
}

