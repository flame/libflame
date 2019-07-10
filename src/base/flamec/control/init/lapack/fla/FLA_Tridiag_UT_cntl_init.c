/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_tridiagut_t*    fla_tridiagut_cntl_fused = NULL;
fla_tridiagut_t*    fla_tridiagut_cntl_nofus = NULL;
fla_tridiagut_t*    fla_tridiagut_cntl_plain = NULL;

fla_blocksize_t*    fla_tridiagut_bsize_leaf = NULL;

void FLA_Tridiag_UT_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_tridiagut_bsize_leaf = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_tridiagut_bsize_leaf, FLA_TRIDIAG_INNER_TO_OUTER_B_RATIO );

	// Create a control tree that uses fused subproblems.
	fla_tridiagut_cntl_fused = FLA_Cntl_tridiagut_obj_create( FLA_FLAT, 
	                                                          FLA_BLK_FUS_VARIANT3,
	                                                          fla_tridiagut_bsize_leaf );

	// Create a control tree that does not used any fusing.
	fla_tridiagut_cntl_nofus = FLA_Cntl_tridiagut_obj_create( FLA_FLAT, 
	                                                          FLA_BLOCKED_VARIANT3,
	                                                          fla_tridiagut_bsize_leaf );

	// Create a control tree that does not used any re-arrangement in house.
	fla_tridiagut_cntl_plain = FLA_Cntl_tridiagut_obj_create( FLA_FLAT, 
	                                                          FLA_BLOCKED_VARIANT1,
	                                                          fla_tridiagut_bsize_leaf );
}

void FLA_Tridiag_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_tridiagut_cntl_fused );
	FLA_Cntl_obj_free( fla_tridiagut_cntl_nofus );
	FLA_Cntl_obj_free( fla_tridiagut_cntl_plain );

	FLA_Blocksize_free( fla_tridiagut_bsize_leaf );
}

