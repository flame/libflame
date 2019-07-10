/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scalr_t* flash_scalr_cntl;
extern fla_gemm_t*  flash_gemm_cntl_pb_bb;

fla_syrk_t*         flash_syrk_cntl_blas = NULL;
fla_syrk_t*         flash_syrk_cntl_ip = NULL;
fla_syrk_t*         flash_syrk_cntl_op = NULL;
fla_syrk_t*         flash_syrk_cntl_mm = NULL;
fla_blocksize_t*    flash_syrk_bsize = NULL;

void FLASH_Syrk_cntl_init()
{
	// Set syrk blocksize for hierarchical storage.
	flash_syrk_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_syrk_cntl_blas  = FLA_Cntl_syrk_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A * A' forms an inner panel product.
	flash_syrk_cntl_ip    = FLA_Cntl_syrk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_syrk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_syrk_cntl_blas,
	                                                  NULL );

	// Create a control tree that assumes A * A' forms an outer panel product.
	flash_syrk_cntl_op    = FLA_Cntl_syrk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT2,
	                                                  flash_syrk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_syrk_cntl_blas,
	                                                  flash_gemm_cntl_pb_bb );

	// Create a control tree that assumes A is large.
	flash_syrk_cntl_mm    = FLA_Cntl_syrk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_syrk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_syrk_cntl_op,
	                                                  NULL );
}

void FLASH_Syrk_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_syrk_cntl_blas );

	FLA_Cntl_obj_free( flash_syrk_cntl_ip );
	FLA_Cntl_obj_free( flash_syrk_cntl_op );
	FLA_Cntl_obj_free( flash_syrk_cntl_mm );

	FLA_Blocksize_free( flash_syrk_bsize );
}

