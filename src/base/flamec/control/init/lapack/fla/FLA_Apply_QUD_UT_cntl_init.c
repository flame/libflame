/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern TLS_CLASS_SPEC fla_trsm_t*  fla_trsm_cntl_blas;
extern TLS_CLASS_SPEC fla_copyt_t* fla_copyt_cntl_blas;
extern TLS_CLASS_SPEC fla_axpyt_t* fla_axpyt_cntl_blas;

TLS_CLASS_SPEC fla_apqudut_t*      fla_apqudut_cntl_leaf = NULL;
TLS_CLASS_SPEC fla_blocksize_t*    fla_apqudut_var1_bsize = NULL;

void FLA_Apply_QUD_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_apqudut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_apqudut_var1_bsize, FLA_UDDATE_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	fla_apqudut_cntl_leaf = FLA_Cntl_apqudut_obj_create( FLA_FLAT,
	                                                     FLA_BLOCKED_VARIANT1, 
	                                                     fla_apqudut_var1_bsize,
	                                                     NULL,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_trsm_cntl_blas,
	                                                     fla_copyt_cntl_blas,
	                                                     fla_axpyt_cntl_blas );
}

void FLA_Apply_QUD_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_apqudut_cntl_leaf );

	FLA_Blocksize_free( fla_apqudut_var1_bsize );
}

