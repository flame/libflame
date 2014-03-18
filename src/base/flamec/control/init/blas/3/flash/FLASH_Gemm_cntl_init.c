
#include "FLAME.h"

extern fla_scal_t* flash_scal_cntl;

fla_gemm_t*      flash_gemm_cntl_blas;
fla_gemm_t*      flash_gemm_cntl_mm_mp;
fla_gemm_t*      flash_gemm_cntl_mm_pm;
fla_gemm_t*      flash_gemm_cntl_mm_op;
fla_gemm_t*      flash_gemm_cntl_mp_pb;
fla_gemm_t*      flash_gemm_cntl_mp_ip;
fla_gemm_t*      flash_gemm_cntl_pm_bp;
fla_gemm_t*      flash_gemm_cntl_pm_ip;
fla_gemm_t*      flash_gemm_cntl_op_bp;
fla_gemm_t*      flash_gemm_cntl_op_pb;
fla_gemm_t*      flash_gemm_cntl_pb_bb;
fla_gemm_t*      flash_gemm_cntl_bp_bb;
fla_gemm_t*      flash_gemm_cntl_ip_bb;

fla_gemm_t*      flash_gemm_cntl_mm;
fla_gemm_t*      flash_gemm_cntl_mp;
fla_gemm_t*      flash_gemm_cntl_pm;
fla_gemm_t*      flash_gemm_cntl_op;
fla_gemm_t*      flash_gemm_cntl_pb;
fla_gemm_t*      flash_gemm_cntl_bp;
fla_gemm_t*      flash_gemm_cntl_ip;

fla_blocksize_t* flash_gemm_bsize;

void FLASH_Gemm_cntl_init()
{
	// Set gemm blocksize for hierarchical storage.
	flash_gemm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree node that executes a gemm subproblem.
	flash_gemm_cntl_blas  = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create control trees for situations where one dimension is large.
	flash_gemm_cntl_pb_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );
	flash_gemm_cntl_bp_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );
	flash_gemm_cntl_ip_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );

	// Create control trees for situations where two dimensions are large.
	flash_gemm_cntl_mp_ip = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_ip_bb );
	flash_gemm_cntl_op_bp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_bp_bb );
	flash_gemm_cntl_pm_ip = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_ip_bb );
	flash_gemm_cntl_op_pb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pb_bb );
	flash_gemm_cntl_mp_pb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pb_bb );
	flash_gemm_cntl_pm_bp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_bp_bb );

	// Create control trees for situations where all dimensions are large.
	flash_gemm_cntl_mm_pm = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pm_ip );
	flash_gemm_cntl_mm_mp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_mp_ip );
	flash_gemm_cntl_mm_op = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_op_bp );

	// Alias select control trees for convenience, when the caller doesn't
	// care (as he usually doesn't when partitioning hierarchical matricies)
	// which order the matrix is partitioned into blocks
	flash_gemm_cntl_mm = flash_gemm_cntl_mm_op;
	flash_gemm_cntl_mp = flash_gemm_cntl_mp_pb;
	flash_gemm_cntl_pm = flash_gemm_cntl_pm_bp;
	flash_gemm_cntl_op = flash_gemm_cntl_op_pb;
	flash_gemm_cntl_pb = flash_gemm_cntl_pb_bb;
	flash_gemm_cntl_bp = flash_gemm_cntl_bp_bb;
	flash_gemm_cntl_ip = flash_gemm_cntl_ip_bb;
	
}

void FLASH_Gemm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_gemm_cntl_blas );

	FLA_Cntl_obj_free( flash_gemm_cntl_pb_bb );
	FLA_Cntl_obj_free( flash_gemm_cntl_bp_bb );
	FLA_Cntl_obj_free( flash_gemm_cntl_ip_bb );

	FLA_Cntl_obj_free( flash_gemm_cntl_mp_ip );
	FLA_Cntl_obj_free( flash_gemm_cntl_op_bp );
	FLA_Cntl_obj_free( flash_gemm_cntl_pm_ip );
	FLA_Cntl_obj_free( flash_gemm_cntl_op_pb );
	FLA_Cntl_obj_free( flash_gemm_cntl_mp_pb );
	FLA_Cntl_obj_free( flash_gemm_cntl_pm_bp );

	FLA_Cntl_obj_free( flash_gemm_cntl_mm_pm );
	FLA_Cntl_obj_free( flash_gemm_cntl_mm_mp );
	FLA_Cntl_obj_free( flash_gemm_cntl_mm_op );

	FLA_Blocksize_free( flash_gemm_bsize );
}

