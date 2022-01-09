/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_gemm_t*  flash_gemm_cntl_op_bp;
extern TLS_CLASS_SPEC fla_trsm_t*  flash_trsm_cntl_bp;
extern TLS_CLASS_SPEC fla_appiv_t* flash_appiv_cntl_bp;

TLS_CLASS_SPEC fla_lu_t*           flash_lu_piv_cntl_leaf = NULL;
TLS_CLASS_SPEC fla_lu_t*           flash_lu_piv_cntl = NULL;
TLS_CLASS_SPEC fla_blocksize_t*    flash_lu_piv_bsize = NULL;

void FLASH_LU_piv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_lu_piv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_lu_piv_cntl_leaf   = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is large.
	flash_lu_piv_cntl        = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   flash_lu_piv_bsize,
	                                                   flash_lu_piv_cntl_leaf,
	                                                   flash_gemm_cntl_op_bp,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_trsm_cntl_bp,
	                                                   flash_trsm_cntl_bp,
	                                                   flash_appiv_cntl_bp,
	                                                   NULL );
}

void FLASH_LU_piv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_lu_piv_cntl_leaf );
	FLA_Cntl_obj_free( flash_lu_piv_cntl );

	FLA_Blocksize_free( flash_lu_piv_bsize );
}

