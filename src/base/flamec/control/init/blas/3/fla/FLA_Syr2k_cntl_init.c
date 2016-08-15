/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scalr_t* fla_scalr_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;

fla_syr2k_t*        fla_syr2k_cntl_blas = NULL;
fla_syr2k_t*        fla_syr2k_cntl_ip = NULL;
fla_syr2k_t*        fla_syr2k_cntl_op = NULL;
fla_syr2k_t*        fla_syr2k_cntl_mm = NULL;
fla_blocksize_t*    fla_syr2k_var3_bsize = NULL;
fla_blocksize_t*    fla_syr2k_var9_bsize = NULL;

void FLA_Syr2k_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_syr2k_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_syr2k_var9_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A and B are b x b blocks.
	fla_syr2k_cntl_blas  = FLA_Cntl_syr2k_obj_create( FLA_FLAT,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A and B form an inner panel product.
	fla_syr2k_cntl_ip    = FLA_Cntl_syr2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  fla_syr2k_var9_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_syr2k_cntl_blas,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A and B form an outer panel product.
	fla_syr2k_cntl_op    = FLA_Cntl_syr2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  fla_syr2k_var3_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_syr2k_cntl_blas,
	                                                  fla_gemm_cntl_blas,
	                                                  fla_gemm_cntl_blas );

	// Create a control tree that assumes A and B are both large.
	fla_syr2k_cntl_mm    = FLA_Cntl_syr2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  fla_syr2k_var9_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_syr2k_cntl_op,
	                                                  NULL,
	                                                  NULL );
}

void FLA_Syr2k_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_syr2k_cntl_blas );

	FLA_Cntl_obj_free( fla_syr2k_cntl_ip );
	FLA_Cntl_obj_free( fla_syr2k_cntl_op );
	FLA_Cntl_obj_free( fla_syr2k_cntl_mm );

	FLA_Blocksize_free( fla_syr2k_var3_bsize );
	FLA_Blocksize_free( fla_syr2k_var9_bsize );
}

