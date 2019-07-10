/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl;
extern fla_gemm_t* flash_gemm_cntl_op_bp;

fla_trsm_t*        flash_trsm_cntl_blas = NULL;
fla_trsm_t*        flash_trsm_cntl_bp = NULL;
fla_trsm_t*        flash_trsm_cntl_mp = NULL;
fla_trsm_t*        flash_trsm_cntl_mm = NULL;
fla_blocksize_t*   flash_trsm_bsize = NULL;

void FLASH_Trsm_cntl_init()
{
	// Set trsm blocksize for hierarchical storage.
	flash_trsm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_trsm_cntl_blas  = FLA_Cntl_trsm_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	flash_trsm_cntl_bp    = FLA_Cntl_trsm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_trsm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trsm_cntl_blas,
	                                                  NULL );

	// Create a control tree that assumes A is large and B is a panel.
	flash_trsm_cntl_mp    = FLA_Cntl_trsm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT2,
	                                                  flash_trsm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trsm_cntl_blas,
	                                                  flash_gemm_cntl_op_bp );

	// Create a control tree that assumes A and B are both large.
	flash_trsm_cntl_mm    = FLA_Cntl_trsm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_trsm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trsm_cntl_mp,
	                                                  NULL );
}

void FLASH_Trsm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_trsm_cntl_blas );

	FLA_Cntl_obj_free( flash_trsm_cntl_bp );
	FLA_Cntl_obj_free( flash_trsm_cntl_mp );
	FLA_Cntl_obj_free( flash_trsm_cntl_mm );

	FLA_Blocksize_free( flash_trsm_bsize );
}

