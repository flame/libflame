/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl;
extern fla_gemm_t* flash_gemm_cntl_op_bp;
extern fla_gemm_t* flash_gemm_cntl_mm_pm;
extern fla_gemm_t* flash_gemm_cntl_mm_mp;

fla_symm_t*        flash_symm_cntl_blas = NULL;
fla_symm_t*        flash_symm_cntl_bp = NULL;
fla_symm_t*        flash_symm_cntl_mp = NULL;
fla_symm_t*        flash_symm_cntl_mm = NULL;
fla_blocksize_t*   flash_symm_bsize = NULL;

void FLASH_Symm_cntl_init()
{
	// Set symm blocksize for hierarchical storage.
	flash_symm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_symm_cntl_blas  = FLA_Cntl_symm_obj_create( FLA_HIER, 
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	flash_symm_cntl_bp    = FLA_Cntl_symm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_symm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_symm_cntl_blas,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is large and B is a panel.
	flash_symm_cntl_mp    = FLA_Cntl_symm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_symm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_symm_cntl_blas,
	                                                  flash_gemm_cntl_op_bp,
	                                                  flash_gemm_cntl_mm_mp );

	// Create a control tree that assumes A and B are both large.
	flash_symm_cntl_mm    = FLA_Cntl_symm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_symm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_symm_cntl_mp,
	                                                  NULL,
	                                                  NULL );
}

void FLASH_Symm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_symm_cntl_blas );

	FLA_Cntl_obj_free( flash_symm_cntl_bp );
	FLA_Cntl_obj_free( flash_symm_cntl_mp );
	FLA_Cntl_obj_free( flash_symm_cntl_mm );

	FLA_Blocksize_free( flash_symm_bsize );
}

