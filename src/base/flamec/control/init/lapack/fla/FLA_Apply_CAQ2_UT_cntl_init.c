/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t* fla_gemm_cntl_blas;
extern fla_trmm_t* fla_trmm_cntl_blas;
extern fla_trsm_t* fla_trsm_cntl_blas;
extern fla_axpy_t* fla_axpy_cntl_blas;
extern fla_copy_t* fla_copy_cntl_blas;

fla_apcaq2ut_t*     fla_apcaq2ut_cntl_leaf = NULL;
fla_blocksize_t*    fla_apcaq2ut_var1_bsize = NULL;

void FLA_Apply_CAQ2_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_apcaq2ut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_apcaq2ut_var1_bsize, FLA_CAQR_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	fla_apcaq2ut_cntl_leaf = FLA_Cntl_apcaq2ut_obj_create( FLA_FLAT,
	                                                       FLA_BLOCKED_VARIANT1, 
	                                                       fla_apcaq2ut_var1_bsize,
	                                                       NULL,
	                                                       fla_gemm_cntl_blas,
	                                                       fla_gemm_cntl_blas,
	                                                       fla_trmm_cntl_blas,
	                                                       fla_trmm_cntl_blas,
	                                                       fla_trsm_cntl_blas,
	                                                       fla_axpy_cntl_blas,
	                                                       fla_axpy_cntl_blas,
	                                                       fla_axpy_cntl_blas,
	                                                       fla_copy_cntl_blas );
}

void FLA_Apply_CAQ2_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_apcaq2ut_cntl_leaf );

	FLA_Blocksize_free( fla_apcaq2ut_var1_bsize );
}

