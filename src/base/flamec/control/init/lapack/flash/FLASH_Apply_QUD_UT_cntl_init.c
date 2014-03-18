
#include "FLAME.h"

fla_apqudut_t*   flash_apqudut_cntl_leaf;
fla_apqudut_t*   flash_apqudut_cntl_mid;
fla_apqudut_t*   flash_apqudut_cntl;
fla_blocksize_t* flash_apqudut_var2_bsize;
fla_blocksize_t* flash_apqudut_var3_bsize;

void FLASH_Apply_QUD_UT_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_apqudut_var2_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );
	flash_apqudut_var3_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_apqudut_cntl_leaf = FLA_Cntl_apqudut_obj_create( FLA_HIER,
	                                                       FLA_SUBPROBLEM, 
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL );

	// Create a control tree to invoke variant 2.
	flash_apqudut_cntl_mid  = FLA_Cntl_apqudut_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT2, 
	                                                       flash_apqudut_var2_bsize,
	                                                       flash_apqudut_cntl_leaf,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL );

	// Create a control tree to invoke variant 3.
	flash_apqudut_cntl      = FLA_Cntl_apqudut_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT3, 
	                                                       flash_apqudut_var3_bsize,
	                                                       flash_apqudut_cntl_mid,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL,
	                                                       NULL );
}

void FLASH_Apply_QUD_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_apqudut_cntl_leaf );
	FLA_Cntl_obj_free( flash_apqudut_cntl_mid );
	FLA_Cntl_obj_free( flash_apqudut_cntl );

	FLA_Blocksize_free( flash_apqudut_var2_bsize );
	FLA_Blocksize_free( flash_apqudut_var3_bsize );
}

