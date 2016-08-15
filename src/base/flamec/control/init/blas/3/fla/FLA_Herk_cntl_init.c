/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scalr_t* fla_scalr_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;

fla_herk_t*         fla_herk_cntl_blas = NULL;
fla_herk_t*         fla_herk_cntl_ip = NULL;
fla_herk_t*         fla_herk_cntl_op = NULL;
fla_herk_t*         fla_herk_cntl_mm = NULL;
fla_blocksize_t*    fla_herk_var2_bsize = NULL;
fla_blocksize_t*    fla_herk_var5_bsize = NULL;

void FLA_Herk_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_herk_var2_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_herk_var5_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A is a b x b block.
	fla_herk_cntl_blas  = FLA_Cntl_herk_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create a control tree that assumes A * A' forms an inner panel product.
	fla_herk_cntl_ip    = FLA_Cntl_herk_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT5,
	                                                fla_herk_var5_bsize,
	                                                fla_scalr_cntl_blas,
	                                                fla_herk_cntl_blas,
	                                                NULL );

	// Create a control tree that assumes A * A' forms an outer panel product.
	fla_herk_cntl_op    = FLA_Cntl_herk_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT2,
	                                                fla_herk_var2_bsize,
	                                                fla_scalr_cntl_blas,
	                                                fla_herk_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create a control tree that assumes A is large.
	fla_herk_cntl_mm    = FLA_Cntl_herk_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT5,
	                                                fla_herk_var5_bsize,
	                                                fla_scalr_cntl_blas,
	                                                fla_herk_cntl_op,
	                                                NULL );
}

void FLA_Herk_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_herk_cntl_blas );

	FLA_Cntl_obj_free( fla_herk_cntl_ip );
	FLA_Cntl_obj_free( fla_herk_cntl_op );
	FLA_Cntl_obj_free( fla_herk_cntl_mm );

	FLA_Blocksize_free( fla_herk_var2_bsize );
	FLA_Blocksize_free( fla_herk_var5_bsize );
}

