/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

fla_apqut_t*        fla_apqut_cntl_leaf = NULL;
fla_apqut_t*        fla_apqut_cntl = NULL;
fla_blocksize_t*    fla_apqut_var1_bsize = NULL;
fla_blocksize_t*    fla_apqut_var2_bsize = NULL;

void FLA_Apply_Q_UT_cntl_init()
{
	// Set the outer blocksize to the default value for conventional storage,
	// and the inner blocksize to the same value but scaled down.
	fla_apqut_var2_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_apqut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_apqut_var1_bsize, FLA_QR_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	fla_apqut_cntl_leaf = FLA_Cntl_apqut_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT1, 
	                                                 fla_apqut_var1_bsize,
	                                                 NULL,
	                                                 fla_trmm_cntl_blas,
	                                                 fla_trmm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_copyt_cntl_blas,
	                                                 fla_axpyt_cntl_blas );
/*
	// Create a control tree to invoke variant 2.
	fla_apqut_cntl      = FLA_Cntl_apqut_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT2, 
	                                                 fla_apqut_var2_bsize,
	                                                 fla_apqut_cntl_leaf,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );
*/
}

void FLA_Apply_Q_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_apqut_cntl_leaf );
/*
	FLA_Cntl_obj_free( fla_apqut_cntl );
*/

	FLA_Blocksize_free( fla_apqut_var1_bsize );
	FLA_Blocksize_free( fla_apqut_var2_bsize );
}

