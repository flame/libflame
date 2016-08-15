/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_uddateut_t*  flash_uddateut_cntl_leaf = NULL;
fla_uddateut_t*  flash_uddateut_cntl = NULL;
fla_blocksize_t* flash_uddateut_var2_bsize = NULL;

void FLASH_UDdate_UT_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_uddateut_var2_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_uddateut_cntl_leaf = FLA_Cntl_uddateut_obj_create( FLA_HIER,
	                                                         FLA_SUBPROBLEM, 
	                                                         NULL,
	                                                         NULL,
	                                                         NULL );

	// Create a control tree to invoke variant 2.
	flash_uddateut_cntl    = FLA_Cntl_uddateut_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT2, 
	                                                       flash_uddateut_var2_bsize,
	                                                       flash_uddateut_cntl_leaf,
	                                                       NULL );
}

void FLASH_UDdate_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_uddateut_cntl_leaf );
	FLA_Cntl_obj_free( flash_uddateut_cntl );

	FLA_Blocksize_free( flash_uddateut_var2_bsize );
}

