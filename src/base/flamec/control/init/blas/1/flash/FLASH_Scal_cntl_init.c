/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_scal_t*        flash_scal_cntl_blas = NULL;
fla_scal_t*        flash_scal_cntl_tb;
fla_scal_t*        flash_scal_cntl_lr;
fla_scal_t*        flash_scal_cntl;
fla_blocksize_t*   flash_scal_bsize;

void FLASH_Scal_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_scal_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is small.
	flash_scal_cntl_blas = FLA_Cntl_scal_obj_create( FLA_HIER,
	                                                 FLA_SUBPROBLEM,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree that marches through A vertically.
	flash_scal_cntl_tb   = FLA_Cntl_scal_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT1,
	                                                 flash_scal_bsize,
	                                                 flash_scal_cntl_blas );

	// Create a control tree that marches through A horizontally.
	flash_scal_cntl_lr   = FLA_Cntl_scal_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT3,
	                                                 flash_scal_bsize,
	                                                 flash_scal_cntl_blas );

	// Create a control tree that marches through A horizontally, then
	// vertically.
	flash_scal_cntl      = FLA_Cntl_scal_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT3,
	                                                 flash_scal_bsize,
	                                                 flash_scal_cntl_tb );
}

void FLASH_Scal_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_scal_cntl_blas );

	FLA_Cntl_obj_free( flash_scal_cntl_tb );
	FLA_Cntl_obj_free( flash_scal_cntl_lr );
	FLA_Cntl_obj_free( flash_scal_cntl );

	FLA_Blocksize_free( flash_scal_bsize );
}

