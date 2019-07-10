/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_caqr2ut_t*  flash_caqr2ut_cntl;
extern fla_apcaq2ut_t* flash_apcaq2ut_cntl;

fla_caqrutinc_t*     flash_caqrutinc_cntl = NULL;
fla_blocksize_t*     flash_caqrutinc_var1_bsize = NULL;

void FLASH_CAQR_UT_inc_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_caqrutinc_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_caqrutinc_cntl = FLA_Cntl_caqrutinc_obj_create( FLA_HIER,
	                                                      FLA_BLOCKED_VARIANT1, 
	                                                      flash_caqrutinc_var1_bsize,
	                                                      flash_caqr2ut_cntl,
	                                                      flash_apcaq2ut_cntl );
}

void FLASH_CAQR_UT_inc_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_caqrutinc_cntl );

	FLA_Blocksize_free( flash_caqrutinc_var1_bsize );
}

