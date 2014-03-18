
#include "FLAME.h"

extern fla_gemm_t* flash_gemm_cntl_pm_bp;
extern fla_gemm_t* flash_gemm_cntl_ip_bb;

fla_sylv_t*        flash_sylv_cntl_leaf;
fla_sylv_t*        flash_sylv_cntl_mb;
fla_sylv_t*        flash_sylv_cntl;
fla_blocksize_t*   flash_sylv_bsize;

void FLASH_Sylv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_sylv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_sylv_cntl_leaf   = FLA_Cntl_sylv_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is a matrix and B is a block.
	flash_sylv_cntl_mb     = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT17,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_leaf,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_ip_bb,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is a matrix and B is a matrix.
	flash_sylv_cntl        = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT15,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_mb,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_pm_bp,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );
}

void FLASH_Sylv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_sylv_cntl_leaf );
	FLA_Cntl_obj_free( flash_sylv_cntl_mb );
	FLA_Cntl_obj_free( flash_sylv_cntl );

	FLA_Blocksize_free( flash_sylv_bsize );
}

