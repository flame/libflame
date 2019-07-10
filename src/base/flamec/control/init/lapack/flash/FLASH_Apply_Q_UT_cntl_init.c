/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_trmm_t*  flash_trmm_cntl_bp;
extern fla_trmm_t*  flash_trmm_cntl_bp;
extern fla_gemm_t*  flash_gemm_cntl_pm;
extern fla_gemm_t*  flash_gemm_cntl_op;
extern fla_trsm_t*  flash_trsm_cntl_bp;
extern fla_copyt_t*  flash_copyt_cntl;
extern fla_axpyt_t*  flash_axpyt_cntl;

fla_apqut_t*        flash_apqut_cntl_leaf = NULL;
fla_apqut_t*        flash_apqut_cntl = NULL;
fla_apqut_t*        flash_apqut_cntl_blas = NULL;
fla_blocksize_t*    flash_apqut_var1_bsize = NULL;
fla_blocksize_t*    flash_apqut_var2_bsize = NULL;

void FLASH_Apply_Q_UT_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_apqut_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );
	flash_apqut_var2_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to dereference block operands and perform
	// flat subproblem.
	flash_apqut_cntl_leaf = FLA_Cntl_apqut_obj_create( FLA_HIER,
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

	// Create a control tree to invoke variant 2 to further partition blocks.
	flash_apqut_cntl    = FLA_Cntl_apqut_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT2, 
	                                                 flash_apqut_var2_bsize,
	                                                 flash_apqut_cntl_leaf,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 3, using hierarchical level-3
	// BLAS control trees.
	flash_apqut_cntl_blas = FLA_Cntl_apqut_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3, 
	                                                   flash_apqut_var1_bsize,
	                                                   NULL,
	                                                   flash_trmm_cntl_bp,
	                                                   flash_trmm_cntl_bp,
	                                                   flash_gemm_cntl_pm,
	                                                   flash_gemm_cntl_op,
	                                                   flash_trsm_cntl_bp,
	                                                   flash_copyt_cntl,
	                                                   flash_axpyt_cntl );
}

void FLASH_Apply_Q_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_apqut_cntl_leaf );
	FLA_Cntl_obj_free( flash_apqut_cntl );
	FLA_Cntl_obj_free( flash_apqut_cntl_blas );

	FLA_Blocksize_free( flash_apqut_var1_bsize );
	FLA_Blocksize_free( flash_apqut_var2_bsize );
}

