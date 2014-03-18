
#include "FLAME.h"

extern fla_gemm_t* flash_gemm_cntl_op_bp;
extern fla_trsm_t* flash_trsm_cntl_bp;

fla_trinv_t*       flash_trinv_cntl_leaf;
fla_trinv_t*       flash_trinv_cntl;
fla_blocksize_t*   flash_trinv_bsize;

void FLASH_Trinv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_trinv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_trinv_cntl_leaf   = FLA_Cntl_trinv_obj_create( FLA_HIER,
	                                                     FLA_SUBPROBLEM,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL,
	                                                     NULL );

	// Create a control tree that assumes A is large.
	flash_trinv_cntl        = FLA_Cntl_trinv_obj_create( FLA_HIER,
	                                                     FLA_BLOCKED_VARIANT3, 
	                                                     flash_trinv_bsize,
	                                                     flash_trinv_cntl_leaf,
	                                                     NULL,
	                                                     flash_trsm_cntl_bp,
	                                                     flash_trsm_cntl_bp,
	                                                     flash_gemm_cntl_op_bp );
}

void FLASH_Trinv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_trinv_cntl_leaf );
	FLA_Cntl_obj_free( flash_trinv_cntl );

	FLA_Blocksize_free( flash_trinv_bsize );
}

