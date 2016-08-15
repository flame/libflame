/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_axpyt_t*       flash_axpyt_cntl_blas = NULL;
fla_axpyt_t*       flash_axpyt_cntl_tb;
fla_axpyt_t*       flash_axpyt_cntl_lr;
fla_axpyt_t*       flash_axpyt_cntl;
fla_blocksize_t*   flash_axpyt_bsize;

void FLASH_Axpyt_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_axpyt_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are small.
	flash_axpyt_cntl_blas = FLA_Cntl_axpyt_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that marches through A and B vertically.
	flash_axpyt_cntl_tb   = FLA_Cntl_axpyt_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   flash_axpyt_bsize,
	                                                   flash_axpyt_cntl_blas );

	// Create a control tree that marches through A and B horizontally.
	flash_axpyt_cntl_lr   = FLA_Cntl_axpyt_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   flash_axpyt_bsize,
	                                                   flash_axpyt_cntl_blas );

	// Create a control tree that marches through A and B horizontally, then
	// vertically.
	flash_axpyt_cntl      = FLA_Cntl_axpyt_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   flash_axpyt_bsize,
	                                                   flash_axpyt_cntl_tb );
}

void FLASH_Axpyt_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_axpyt_cntl_blas );

	FLA_Cntl_obj_free( flash_axpyt_cntl_tb );
	FLA_Cntl_obj_free( flash_axpyt_cntl_lr );
	FLA_Cntl_obj_free( flash_axpyt_cntl );

	FLA_Blocksize_free( flash_axpyt_bsize );
}

