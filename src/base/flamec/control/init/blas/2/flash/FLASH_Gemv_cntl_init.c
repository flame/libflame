/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl;

fla_gemv_t*      flash_gemv_cntl_blas = NULL;
fla_gemv_t*      flash_gemv_cntl_fm_rp;
fla_gemv_t*      flash_gemv_cntl_fm_cp;
fla_gemv_t*      flash_gemv_cntl_rp_bv;
fla_gemv_t*      flash_gemv_cntl_cp_bv;
fla_blocksize_t* flash_gemv_bsize;

void FLASH_Gemv_cntl_init()
{
	// Set gemv blocksize for hierarchical storage.
	flash_gemv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree node that executes a gemv subproblem.
	flash_gemv_cntl_blas   = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create control trees for situations where one dimension is large.
	flash_gemv_cntl_cp_bv = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_blas );
	flash_gemv_cntl_rp_bv = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_blas );

	// Create control trees for situations where both dimensions are large.
	flash_gemv_cntl_fm_rp = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_rp_bv );
	flash_gemv_cntl_fm_cp = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_cp_bv );
}

void FLASH_Gemv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_gemv_cntl_blas );

	FLA_Cntl_obj_free( flash_gemv_cntl_cp_bv );
	FLA_Cntl_obj_free( flash_gemv_cntl_rp_bv );

	FLA_Cntl_obj_free( flash_gemv_cntl_fm_rp );
	FLA_Cntl_obj_free( flash_gemv_cntl_fm_cp );

	FLA_Blocksize_free( flash_gemv_bsize );
}

