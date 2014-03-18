
#include "FLAME.h"

extern fla_apqut_t*  flash_apqut_cntl_blas;

fla_qrut_t*          flash_qrut_cntl_leaf;
fla_qrut_t*          flash_qrut_cntl;

fla_blocksize_t*     flash_qrut_var3_bsize;

void FLASH_QR_UT_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_qrut_var3_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to compute the subproblem.
	flash_qrut_cntl_leaf = FLA_Cntl_qrut_obj_create( FLA_HIER,
	                                                 FLA_SUBPROBLEM, 
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 3.
	flash_qrut_cntl      = FLA_Cntl_qrut_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT3,
	                                                 flash_qrut_var3_bsize,
	                                                 flash_qrut_cntl_leaf,
	                                                 flash_apqut_cntl_blas );
}

void FLASH_QR_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_qrut_cntl_leaf );
	FLA_Cntl_obj_free( flash_qrut_cntl );

	FLA_Blocksize_free( flash_qrut_var3_bsize );
}

