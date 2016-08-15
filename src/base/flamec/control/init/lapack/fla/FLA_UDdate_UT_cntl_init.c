/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_apqudut_t* fla_apqudut_cntl_leaf;

fla_uddateut_t*       fla_uddateut_cntl_unb = NULL;
fla_uddateut_t*       fla_uddateut_cntl_leaf = NULL;
fla_blocksize_t*      fla_uddateut_var1_bsize = NULL;

void FLA_UDdate_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_uddateut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_uddateut_var1_bsize, FLA_UDDATE_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke unblocked variant 1.
	fla_uddateut_cntl_unb = FLA_Cntl_uddateut_obj_create( FLA_FLAT,
	                                                      FLA_UNB_OPT_VARIANT1, 
	                                                      NULL,
	                                                      NULL,
	                                                      NULL );

	// Create a control tree for small-to-medium sequential problems and
	// as the means to compute on FLASH blocks.
	fla_uddateut_cntl_leaf = FLA_Cntl_uddateut_obj_create( FLA_FLAT,
	                                                       FLA_BLOCKED_VARIANT1, 
	                                                       fla_uddateut_var1_bsize,
	                                                       fla_uddateut_cntl_unb,
	                                                       fla_apqudut_cntl_leaf );

}

void FLA_UDdate_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_uddateut_cntl_unb );
	FLA_Cntl_obj_free( fla_uddateut_cntl_leaf );

	FLA_Blocksize_free( fla_uddateut_var1_bsize );
}

