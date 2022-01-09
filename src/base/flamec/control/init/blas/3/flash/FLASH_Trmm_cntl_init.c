/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_scal_t* flash_scal_cntl;
extern TLS_CLASS_SPEC fla_gemm_t* flash_gemm_cntl_op_bp;

TLS_CLASS_SPEC fla_trmm_t*        flash_trmm_cntl_blas = NULL;
TLS_CLASS_SPEC fla_trmm_t*        flash_trmm_cntl_bp = NULL;
TLS_CLASS_SPEC fla_trmm_t*        flash_trmm_cntl_mp = NULL;
TLS_CLASS_SPEC fla_trmm_t*        flash_trmm_cntl_mm = NULL;
TLS_CLASS_SPEC fla_blocksize_t*   flash_trmm_bsize = NULL;

void FLASH_Trmm_cntl_init()
{
	// Set trmm blocksize for hierarchical storage.
	flash_trmm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_trmm_cntl_blas  = FLA_Cntl_trmm_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	flash_trmm_cntl_bp    = FLA_Cntl_trmm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_trmm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trmm_cntl_blas,
	                                                  NULL );

	// Create a control tree that assumes A is large and B is a panel.
	flash_trmm_cntl_mp    = FLA_Cntl_trmm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT2,
	                                                  flash_trmm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trmm_cntl_blas,
	                                                  flash_gemm_cntl_op_bp );

	// Create a control tree that assumes A and B are both large.
	flash_trmm_cntl_mm    = FLA_Cntl_trmm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_trmm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_trmm_cntl_mp,
	                                                  NULL );
}

void FLASH_Trmm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_trmm_cntl_blas );

	FLA_Cntl_obj_free( flash_trmm_cntl_bp );
	FLA_Cntl_obj_free( flash_trmm_cntl_mp );
	FLA_Cntl_obj_free( flash_trmm_cntl_mm );

	FLA_Blocksize_free( flash_trmm_bsize );
}

