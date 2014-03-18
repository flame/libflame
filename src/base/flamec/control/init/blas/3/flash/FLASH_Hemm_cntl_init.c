
#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl;
extern fla_gemm_t* flash_gemm_cntl_op_bp;
extern fla_gemm_t* flash_gemm_cntl_mm_pm;
extern fla_gemm_t* flash_gemm_cntl_mm_mp;

fla_hemm_t*        flash_hemm_cntl_blas;
fla_hemm_t*        flash_hemm_cntl_bp;
fla_hemm_t*        flash_hemm_cntl_mp;
fla_hemm_t*        flash_hemm_cntl_mm;
fla_blocksize_t*   flash_hemm_bsize;

void FLASH_Hemm_cntl_init()
{
	// Set hemm blocksize for hierarchical storage.
	flash_hemm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_hemm_cntl_blas  = FLA_Cntl_hemm_obj_create( FLA_HIER, 
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	flash_hemm_cntl_bp    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_blas,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is large and B is a panel.
	flash_hemm_cntl_mp    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_blas,
	                                                  flash_gemm_cntl_op_bp,
	                                                  flash_gemm_cntl_mm_mp );

	// Create a control tree that assumes A and B are both large.
	flash_hemm_cntl_mm    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_mp,
	                                                  NULL,
	                                                  NULL );
}

void FLASH_Hemm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_hemm_cntl_blas );

	FLA_Cntl_obj_free( flash_hemm_cntl_bp );
	FLA_Cntl_obj_free( flash_hemm_cntl_mp );
	FLA_Cntl_obj_free( flash_hemm_cntl_mm );

	FLA_Blocksize_free( flash_hemm_bsize );
}

