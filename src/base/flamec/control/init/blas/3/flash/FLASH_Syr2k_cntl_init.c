
#include "FLAME.h"

extern fla_scalr_t* flash_scalr_cntl;
extern fla_gemm_t*  flash_gemm_cntl_pb_bb;

fla_syr2k_t*        flash_syr2k_cntl_blas;
fla_syr2k_t*        flash_syr2k_cntl_ip;
fla_syr2k_t*        flash_syr2k_cntl_op;
fla_syr2k_t*        flash_syr2k_cntl_mm;
fla_blocksize_t*    flash_syr2k_bsize;

void FLASH_Syr2k_cntl_init()
{
	// Set syr2k blocksize for hierarchical storage.
	flash_syr2k_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_syr2k_cntl_blas  = FLA_Cntl_syr2k_obj_create( FLA_HIER,
	                                                    FLA_SUBPROBLEM,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL );

	// Create a control tree that assumes A and B form an inner panel product.
	flash_syr2k_cntl_ip    = FLA_Cntl_syr2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    flash_syr2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_syr2k_cntl_blas,
	                                                    NULL,
	                                                    NULL );

	// Create a control tree that assumes A and B form an outer panel product.
	flash_syr2k_cntl_op    = FLA_Cntl_syr2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT4,
	                                                    flash_syr2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_syr2k_cntl_blas,
	                                                    flash_gemm_cntl_pb_bb,
	                                                    flash_gemm_cntl_pb_bb );

	// Create a control tree that assumes A and B are both large.
	flash_syr2k_cntl_mm    = FLA_Cntl_syr2k_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT9,
	                                                    flash_syr2k_bsize,
	                                                    flash_scalr_cntl,
	                                                    flash_syr2k_cntl_op,
	                                                    NULL,
	                                                    NULL );
}

void FLASH_Syr2k_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_syr2k_cntl_blas );

	FLA_Cntl_obj_free( flash_syr2k_cntl_ip );
	FLA_Cntl_obj_free( flash_syr2k_cntl_op );
	FLA_Cntl_obj_free( flash_syr2k_cntl_mm );

	FLA_Blocksize_free( flash_syr2k_bsize );
}

