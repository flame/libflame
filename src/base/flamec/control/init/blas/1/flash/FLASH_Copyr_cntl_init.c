
#include "FLAME.h"

extern fla_copy_t* flash_copy_cntl_tb;

fla_copyr_t*       flash_copyr_cntl_blas;
fla_copyr_t*       flash_copyr_cntl;
fla_blocksize_t*   flash_copyr_bsize;

void FLASH_Copyr_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_copyr_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are small.
	flash_copyr_cntl_blas = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that marches through A and B from TL to BR.
	flash_copyr_cntl      = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   flash_copyr_bsize,
	                                                   flash_copyr_cntl_blas,
	                                                   flash_copy_cntl_tb );
}

void FLASH_Copyr_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_copyr_cntl_blas );

	FLA_Cntl_obj_free( flash_copyr_cntl );

	FLA_Blocksize_free( flash_copyr_bsize );
}

