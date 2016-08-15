/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_bidiagut_t*    fla_bidiagut_cntl_fused = NULL;
fla_bidiagut_t*    fla_bidiagut_cntl_nofus = NULL;
fla_bidiagut_t*    fla_bidiagut_cntl_plain = NULL;

fla_blocksize_t*   fla_bidiagut_bsize_leaf = NULL;

void FLA_Bidiag_UT_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_bidiagut_bsize_leaf = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_bidiagut_bsize_leaf, FLA_BIDIAG_INNER_TO_OUTER_B_RATIO );

	// Create a control tree that uses fused subproblems.
	fla_bidiagut_cntl_fused = FLA_Cntl_bidiagut_obj_create( FLA_FLAT, 
	                                                        FLA_BLK_FUS_VARIANT4,
	                                                        fla_bidiagut_bsize_leaf );

	// Create a control tree that does not used any fusing.
	fla_bidiagut_cntl_nofus = FLA_Cntl_bidiagut_obj_create( FLA_FLAT, 
	                                                        FLA_BLOCKED_VARIANT4,
	                                                        fla_bidiagut_bsize_leaf );

	// Create a control tree that reflects the basic algorithm.
	fla_bidiagut_cntl_plain = FLA_Cntl_bidiagut_obj_create( FLA_FLAT, 
	                                                        FLA_BLOCKED_VARIANT1,
	                                                        fla_bidiagut_bsize_leaf );
}

void FLA_Bidiag_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_bidiagut_cntl_fused );
	FLA_Cntl_obj_free( fla_bidiagut_cntl_nofus );
	FLA_Cntl_obj_free( fla_bidiagut_cntl_plain );

	FLA_Blocksize_free( fla_bidiagut_bsize_leaf );
}

