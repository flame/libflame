/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLASH_Her2k_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	// Set her2k blocksize for hierarchical storage.
	FLA_cntl_flash_init_i->flash_her2k_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	FLA_cntl_flash_init_i->flash_her2k_cntl_blas  = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_SUBPROBLEM,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL );
 
	// Create a control tree that assumes A and B form an inner panel product.
	FLA_cntl_flash_init_i->flash_her2k_cntl_ip    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    FLA_cntl_flash_init_i->flash_her2k_bsize,
	                                                    FLA_cntl_flash_init_i->flash_scalr_cntl,
	                                                    FLA_cntl_flash_init_i->flash_her2k_cntl_blas,
	                                                    NULL,
	                                                    NULL );

	// Create a control tree that assumes A and B form an outer panel product.
	FLA_cntl_flash_init_i->flash_her2k_cntl_op    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT4,
	                                                    FLA_cntl_flash_init_i->flash_her2k_bsize,
	                                                    FLA_cntl_flash_init_i->flash_scalr_cntl,
	                                                    FLA_cntl_flash_init_i->flash_her2k_cntl_blas,
	                                                    FLA_cntl_flash_init_i->flash_gemm_cntl_pb_bb,
	                                                    FLA_cntl_flash_init_i->flash_gemm_cntl_pb_bb );

	// Create a control tree that assumes A and B are both large.
	FLA_cntl_flash_init_i->flash_her2k_cntl_mm    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    FLA_cntl_flash_init_i->flash_her2k_bsize,
	                                                    FLA_cntl_flash_init_i->flash_scalr_cntl,
	                                                    FLA_cntl_flash_init_i->flash_her2k_cntl_op,
	                                                    NULL,
	                                                    NULL );
}

void FLASH_Her2k_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_her2k_cntl_blas );

	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_her2k_cntl_ip );
	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_her2k_cntl_op );
	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_her2k_cntl_mm );

	FLA_Blocksize_free( FLA_cntl_flash_init_i->flash_her2k_bsize );
}

#endif

extern fla_scalr_t* flash_scalr_cntl;
extern fla_gemm_t*  flash_gemm_cntl_pb_bb;

fla_her2k_t*        flash_her2k_cntl_blas = NULL;
fla_her2k_t*        flash_her2k_cntl_ip = NULL;
fla_her2k_t*        flash_her2k_cntl_op = NULL;
fla_her2k_t*        flash_her2k_cntl_mm = NULL;
fla_blocksize_t*    flash_her2k_bsize = NULL;

void FLASH_Her2k_cntl_init()
{
	// Set her2k blocksize for hierarchical storage.
	flash_her2k_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_her2k_cntl_blas  = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_SUBPROBLEM,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL );
 
	// Create a control tree that assumes A and B form an inner panel product.
	flash_her2k_cntl_ip    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    flash_her2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_her2k_cntl_blas,
	                                                    NULL,
	                                                    NULL );

	// Create a control tree that assumes A and B form an outer panel product.
	flash_her2k_cntl_op    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT4,
	                                                    flash_her2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_her2k_cntl_blas,
	                                                    flash_gemm_cntl_pb_bb,
	                                                    flash_gemm_cntl_pb_bb );

	// Create a control tree that assumes A and B are both large.
	flash_her2k_cntl_mm    = FLA_Cntl_her2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    flash_her2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_her2k_cntl_op,
	                                                    NULL,
	                                                    NULL );
}

void FLASH_Her2k_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_her2k_cntl_blas );

	FLA_Cntl_obj_free( flash_her2k_cntl_ip );
	FLA_Cntl_obj_free( flash_her2k_cntl_op );
	FLA_Cntl_obj_free( flash_her2k_cntl_mm );

	FLA_Blocksize_free( flash_her2k_bsize );
}

