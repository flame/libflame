/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_uddateut_t* flash_uddateut_cntl;
extern fla_apqudut_t*  flash_apqudut_cntl;

fla_uddateutinc_t*     flash_uddateutinc_cntl = NULL;
fla_blocksize_t*       flash_uddateutinc_var1_bsize = NULL;

void FLASH_UDdate_UT_inc_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_uddateutinc_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_uddateutinc_cntl = FLA_Cntl_uddateutinc_obj_create( FLA_HIER,
	                                                          FLA_BLOCKED_VARIANT1, 
	                                                          flash_uddateutinc_var1_bsize,
	                                                          flash_uddateut_cntl,
	                                                          flash_apqudut_cntl );
}

void FLASH_UDdate_UT_inc_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_uddateutinc_cntl );

	FLA_Blocksize_free( flash_uddateutinc_var1_bsize );
}

