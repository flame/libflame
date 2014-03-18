
#include "FLAME.h"

fla_uddateut_t*  flash_uddateut_cntl_leaf;
fla_uddateut_t*  flash_uddateut_cntl;
fla_blocksize_t* flash_uddateut_var2_bsize;

void FLASH_UDdate_UT_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_uddateut_var2_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_uddateut_cntl_leaf = FLA_Cntl_uddateut_obj_create( FLA_HIER,
	                                                         FLA_SUBPROBLEM, 
	                                                         NULL,
	                                                         NULL,
	                                                         NULL );

	// Create a control tree to invoke variant 2.
	flash_uddateut_cntl    = FLA_Cntl_uddateut_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT2, 
	                                                       flash_uddateut_var2_bsize,
	                                                       flash_uddateut_cntl_leaf,
	                                                       NULL );
}

void FLASH_UDdate_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_uddateut_cntl_leaf );
	FLA_Cntl_obj_free( flash_uddateut_cntl );

	FLA_Blocksize_free( flash_uddateut_var2_bsize );
}

