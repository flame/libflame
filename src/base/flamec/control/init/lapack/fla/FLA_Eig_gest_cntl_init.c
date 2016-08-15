/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_axpy_t*  fla_axpy_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;

fla_eig_gest_t*     fla_eig_gest_ix_cntl = NULL;
fla_eig_gest_t*     fla_eig_gest_nx_cntl = NULL;
fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf = NULL;
fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf = NULL;
fla_blocksize_t*    fla_eig_gest_var1_bsize = NULL;

void FLA_Eig_gest_cntl_init()
{
	// Set blocksize with default values for conventional storage.
	fla_eig_gest_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke the unblocked subproblem (inverse cases).
	fla_eig_gest_ix_cntl_leaf = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                          FLA_BLOCKED_EXTERN,
#else
	                                                          FLA_UNB_OPT_VARIANT3,
#endif
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL );

	// Create a control tree to invoke the unblocked subproblem (no inverse cases).
	fla_eig_gest_nx_cntl_leaf = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                          FLA_BLOCKED_EXTERN,
#else
	                                                          FLA_UNB_OPT_VARIANT2,
#endif
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL );

	// Create a control tree for large problems with no extra recursion
	// (inverse cases).
	fla_eig_gest_ix_cntl      = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
	                                                          FLA_BLOCKED_VARIANT4, 
	                                                          fla_eig_gest_var1_bsize,
	                                                          fla_eig_gest_ix_cntl_leaf,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_hemm_cntl_blas,
	                                                          fla_her2k_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trsm_cntl_blas,
	                                                          fla_trsm_cntl_blas );

	// Create a control tree for large problems with no extra recursion
	// (no inverse cases).
	fla_eig_gest_nx_cntl      = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
	                                                          FLA_BLOCKED_VARIANT2, 
	                                                          fla_eig_gest_var1_bsize,
	                                                          fla_eig_gest_nx_cntl_leaf,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_hemm_cntl_blas,
	                                                          fla_her2k_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trsm_cntl_blas,
	                                                          fla_trsm_cntl_blas );
}

void FLA_Eig_gest_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_eig_gest_ix_cntl );
	FLA_Cntl_obj_free( fla_eig_gest_nx_cntl );
	FLA_Cntl_obj_free( fla_eig_gest_ix_cntl_leaf );
	FLA_Cntl_obj_free( fla_eig_gest_nx_cntl_leaf );

	FLA_Blocksize_free( fla_eig_gest_var1_bsize );
}

