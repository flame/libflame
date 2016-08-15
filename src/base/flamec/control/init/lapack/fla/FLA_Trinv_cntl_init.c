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

fla_trinv_t*       fla_trinv_cntl_leaf = NULL;
fla_trinv_t*       fla_trinv_cntl = NULL;
fla_blocksize_t*   fla_trinv_var3_bsize = NULL;

void FLA_Trinv_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_trinv_var3_bsize  = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	//fla_trinv_var3_bsize  = FLA_Blocksize_create( 192, 192, 192, 192 );

	// Create a control tree to invoke LAPACK.
	fla_trinv_cntl_leaf   = FLA_Cntl_trinv_obj_create( FLA_FLAT,
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
	                                                   NULL );

	// Create a control tree to invoke variant 3.
	fla_trinv_cntl        = FLA_Cntl_trinv_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3, 
	                                                   fla_trinv_var3_bsize,
	                                                   fla_trinv_cntl_leaf,
	                                                   fla_trmm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_gemm_cntl_blas );
}

void FLA_Trinv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_trinv_cntl_leaf );
	FLA_Cntl_obj_free( fla_trinv_cntl );

	FLA_Blocksize_free( fla_trinv_var3_bsize );
}

