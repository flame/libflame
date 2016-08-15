/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t* fla_scal_cntl_blas;

fla_gemm_t*      fla_gemm_cntl_blas = NULL;

fla_gemm_t*      fla_gemm_cntl_pb_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_bp_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_ip_bb = NULL;

fla_gemm_t*      fla_gemm_cntl_mp_ip = NULL;
fla_gemm_t*      fla_gemm_cntl_mp_ip_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_op_bp = NULL;
fla_gemm_t*      fla_gemm_cntl_op_bp_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_pm_ip = NULL;
fla_gemm_t*      fla_gemm_cntl_pm_ip_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_op_pb = NULL;
fla_gemm_t*      fla_gemm_cntl_op_pb_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_mp_pb = NULL;
fla_gemm_t*      fla_gemm_cntl_mp_pb_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_pm_bp = NULL;
fla_gemm_t*      fla_gemm_cntl_pm_bp_bb = NULL;

fla_gemm_t*      fla_gemm_cntl_mm_pm = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_pm_ip = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_pm_ip_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_mp = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_mp_ip = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_mp_ip_bb = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_op = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_op_bp = NULL;
fla_gemm_t*      fla_gemm_cntl_mm_op_bp_bb = NULL;

fla_blocksize_t* fla_gemm_var1_bsize = NULL;
fla_blocksize_t* fla_gemm_var3_bsize = NULL;
fla_blocksize_t* fla_gemm_var5_bsize = NULL;

void FLA_Gemm_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_gemm_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_gemm_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_gemm_var5_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree node that executes a gemm subproblem.
	fla_gemm_cntl_blas  = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create control trees for situations where one dimension is large.
	fla_gemm_cntl_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT1,
	                                                fla_gemm_var1_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );
	fla_gemm_cntl_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_gemm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );
	fla_gemm_cntl_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT5,
	                                                fla_gemm_var5_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create control trees for situations where two dimensions are large.
	fla_gemm_cntl_mp_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_mp_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_ip_bb );
	fla_gemm_cntl_op_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_op_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_bp_bb );
	fla_gemm_cntl_pm_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_pm_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_ip_bb );
	fla_gemm_cntl_op_pb    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_op_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_pb_bb );
	fla_gemm_cntl_mp_pb    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_mp_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_pb_bb );
	fla_gemm_cntl_pm_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_pm_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_bp_bb );

	// Create control trees for situations where all dimensions are large.
	fla_gemm_cntl_mm_pm       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_pm_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_pm_ip );
	fla_gemm_cntl_mm_pm_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_pm_ip_bb );
	fla_gemm_cntl_mm_mp       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_mp_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_mp_ip );
	fla_gemm_cntl_mm_mp_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_mp_ip_bb );
	fla_gemm_cntl_mm_op       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_op_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_op_bp );
	fla_gemm_cntl_mm_op_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_op_bp_bb );
}

void FLA_Gemm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_gemm_cntl_blas );

	FLA_Cntl_obj_free( fla_gemm_cntl_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_bp_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_ip_bb );

	FLA_Cntl_obj_free( fla_gemm_cntl_mp_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_bp_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_pb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_pb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_bp_bb );

	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op_bp_bb );

	FLA_Blocksize_free( fla_gemm_var1_bsize );
	FLA_Blocksize_free( fla_gemm_var3_bsize );
	FLA_Blocksize_free( fla_gemm_var5_bsize );
}

