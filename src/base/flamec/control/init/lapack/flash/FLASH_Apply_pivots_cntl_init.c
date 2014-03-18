
#include "FLAME.h"

fla_appiv_t*       flash_appiv_cntl_leaf;
fla_appiv_t*       flash_appiv_cntl_bp;
fla_appiv_t*       flash_appiv_cntl;
fla_blocksize_t*   flash_appiv_bsize;

void FLASH_Apply_pivots_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_appiv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_appiv_cntl_leaf   = FLA_Cntl_appiv_obj_create( FLA_HIER,
	                                                     FLA_SUBPROBLEM,
	                                                     NULL,
	                                                     NULL );

	// Create a control tree that assumes A is large.
	flash_appiv_cntl_bp     = FLA_Cntl_appiv_obj_create( FLA_HIER,
	                                                     FLA_BLOCKED_VARIANT1,
	                                                     flash_appiv_bsize,
	                                                     flash_appiv_cntl_leaf );

	// Create a control tree that assumes A and p are large.
	flash_appiv_cntl        = FLA_Cntl_appiv_obj_create( FLA_HIER,
	                                                     FLA_BLOCKED_VARIANT2,
	                                                     flash_appiv_bsize,
	                                                     flash_appiv_cntl_bp );
}

void FLASH_Apply_pivots_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_appiv_cntl_leaf );
	FLA_Cntl_obj_free( flash_appiv_cntl_bp );
	FLA_Cntl_obj_free( flash_appiv_cntl );

	FLA_Blocksize_free( flash_appiv_bsize );
}

