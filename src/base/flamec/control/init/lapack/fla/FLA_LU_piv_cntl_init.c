/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_appiv_t* fla_appiv_cntl_leaf;

fla_lu_t*           fla_lu_piv_cntl = NULL;
fla_lu_t*           fla_lu_piv_cntl2 = NULL;

fla_lu_t*           fla_lu_piv_cntl_in = NULL;
fla_lu_t*           fla_lu_piv_cntl_leaf = NULL;
fla_blocksize_t*    fla_lu_piv_var5_bsize = NULL;
fla_blocksize_t*    fla_lu_piv_var5_bsize_in = NULL;
double              fla_lu_piv_var5_in_to_ou_bsize_ratio = 0.125;

void FLA_LU_piv_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_lu_piv_var5_bsize    = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_lu_piv_var5_bsize_in = FLA_Blocksize_create_copy( fla_lu_piv_var5_bsize );
	FLA_Blocksize_scale( fla_lu_piv_var5_bsize_in, fla_lu_piv_var5_in_to_ou_bsize_ratio );

	// Create a control tree to invoke LAPACK.
	fla_lu_piv_cntl_leaf   = FLA_Cntl_lu_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                 FLA_BLOCKED_EXTERN,
#else
	                                                 FLA_UNB_OPT_VARIANT4,
#endif
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree for small subproblems.
	fla_lu_piv_cntl_in     = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT5,
	                                                 fla_lu_piv_var5_bsize_in,
	                                                 fla_lu_piv_cntl_leaf,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_appiv_cntl_leaf,
	                                                 fla_appiv_cntl_leaf );

	// Create a control tree for larger problems with one level of recursion.
	fla_lu_piv_cntl2       = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT5,
	                                                 fla_lu_piv_var5_bsize,
	                                                 fla_lu_piv_cntl_in,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_appiv_cntl_leaf,
	                                                 fla_appiv_cntl_leaf );

	// Create a control tree for large problems with no extra recursion.
	fla_lu_piv_cntl        = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT5,
	                                                 fla_lu_piv_var5_bsize,
	                                                 fla_lu_piv_cntl_leaf,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_appiv_cntl_leaf,
	                                                 fla_appiv_cntl_leaf );
}

void FLA_LU_piv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_lu_piv_cntl );
	FLA_Cntl_obj_free( fla_lu_piv_cntl2 );
	FLA_Cntl_obj_free( fla_lu_piv_cntl_in );
	FLA_Cntl_obj_free( fla_lu_piv_cntl_leaf );

	FLA_Blocksize_free( fla_lu_piv_var5_bsize );
	FLA_Blocksize_free( fla_lu_piv_var5_bsize_in );
}

