/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLA_Apply_Q2_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	FLA_cntl_flamec_init_i->fla_apq2ut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale_ts( FLA_cntl_flamec_init_i, FLA_cntl_flamec_init_i->fla_apq2ut_var1_bsize, FLA_QR_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	FLA_cntl_flamec_init_i->fla_apq2ut_cntl_leaf = FLA_Cntl_apq2ut_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1, 
	                                                   FLA_cntl_flamec_init_i->fla_apq2ut_var1_bsize,
	                                                   NULL,
	                                                   FLA_cntl_flamec_init_i->fla_gemm_cntl_blas,
	                                                   FLA_cntl_flamec_init_i->fla_gemm_cntl_blas,
	                                                   FLA_cntl_flamec_init_i->fla_trsm_cntl_blas,
	                                                   FLA_cntl_flamec_init_i->fla_copyt_cntl_blas,
	                                                   FLA_cntl_flamec_init_i->fla_axpyt_cntl_blas );
}

void FLA_Apply_Q2_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flamec_init_i->fla_apq2ut_cntl_leaf );

	FLA_Blocksize_free( FLA_cntl_flamec_init_i->fla_apq2ut_var1_bsize );
}
#endif

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

fla_apq2ut_t*       fla_apq2ut_cntl_leaf = NULL;
fla_blocksize_t*    fla_apq2ut_var1_bsize = NULL;

void FLA_Apply_Q2_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_apq2ut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_apq2ut_var1_bsize, FLA_QR_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	fla_apq2ut_cntl_leaf = FLA_Cntl_apq2ut_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1, 
	                                                   fla_apq2ut_var1_bsize,
	                                                   NULL,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_copyt_cntl_blas,
	                                                   fla_axpyt_cntl_blas );
}

void FLA_Apply_Q2_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_apq2ut_cntl_leaf );

	FLA_Blocksize_free( fla_apq2ut_var1_bsize );
}

