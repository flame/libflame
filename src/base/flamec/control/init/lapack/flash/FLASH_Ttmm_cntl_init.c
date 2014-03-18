
#include "FLAME.h"

extern fla_herk_t* flash_herk_cntl_op;
extern fla_trmm_t* flash_trmm_cntl_bp;

fla_ttmm_t*        flash_ttmm_cntl_leaf;
fla_ttmm_t*        flash_ttmm_cntl;
fla_blocksize_t*   flash_ttmm_bsize;

void FLASH_Ttmm_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_ttmm_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_ttmm_cntl_leaf   = FLA_Cntl_ttmm_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is large.
	flash_ttmm_cntl        = FLA_Cntl_ttmm_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT1, 
	                                                   flash_ttmm_bsize,
	                                                   flash_ttmm_cntl_leaf,
	                                                   flash_herk_cntl_op,
	                                                   flash_trmm_cntl_bp,
	                                                   NULL );
}

void FLASH_Ttmm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_ttmm_cntl_leaf );
	FLA_Cntl_obj_free( flash_ttmm_cntl );

	FLA_Blocksize_free( flash_ttmm_bsize );
}

