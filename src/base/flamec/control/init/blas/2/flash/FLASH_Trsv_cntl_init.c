/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemv_t* flash_gemv_cntl_cp_bv;

fla_trsv_t*        flash_trsv_cntl_blas = NULL;
fla_trsv_t*        flash_trsv_cntl;
fla_blocksize_t*   flash_trsv_bsize;

void FLASH_Trsv_cntl_init()
{
	// Set trsv blocksize for hierarchical storage.
	flash_trsv_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_trsv_cntl_blas  = FLA_Cntl_trsv_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is large.
	flash_trsv_cntl       = FLA_Cntl_trsv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT2,
	                                                  flash_trsv_bsize,
	                                                  flash_trsv_cntl_blas,
	                                                  flash_gemv_cntl_cp_bv );
}

void FLASH_Trsv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_trsv_cntl_blas );

	FLA_Cntl_obj_free( flash_trsv_cntl );

	FLA_Blocksize_free( flash_trsv_bsize );
}

