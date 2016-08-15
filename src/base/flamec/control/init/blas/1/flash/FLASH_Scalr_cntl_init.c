/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl_tb;

fla_scalr_t*       flash_scalr_cntl_blas = NULL;
fla_scalr_t*       flash_scalr_cntl;
fla_blocksize_t*   flash_scalr_bsize;

void FLASH_Scalr_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_scalr_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is small.
	flash_scalr_cntl_blas = FLA_Cntl_scalr_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that computes column panels, top-left to
	// bottom-right.
	flash_scalr_cntl      = FLA_Cntl_scalr_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   flash_scalr_bsize,
	                                                   flash_scalr_cntl_blas,
	                                                   flash_scal_cntl_tb );
}

void FLASH_Scalr_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_scalr_cntl_blas );

	FLA_Cntl_obj_free( flash_scalr_cntl );

	FLA_Blocksize_free( flash_scalr_bsize );
}

