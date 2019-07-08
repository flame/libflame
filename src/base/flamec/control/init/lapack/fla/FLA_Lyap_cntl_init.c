/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLA_Lyap_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	// Set blocksize with default value for conventional storage.
	FLA_cntl_flamec_init_i->fla_lyap_bsize       = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke an unblocked variant.
	FLA_cntl_flamec_init_i->fla_lyap_cntl_leaf   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                                 FLA_UNBLOCKED_VARIANT1,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke a blocked variant.
	FLA_cntl_flamec_init_i->fla_lyap_cntl        = FLA_Cntl_lyap_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT1,
	                                                 FLA_cntl_flamec_init_i->fla_lyap_bsize,
	                                                 FLA_cntl_flamec_init_i->fla_scal_cntl_blas,
	                                                 FLA_cntl_flamec_init_i->fla_lyap_cntl_leaf,
	                                                 FLA_cntl_flamec_init_i->fla_sylv_cntl,
	                                                 FLA_cntl_flamec_init_i->fla_gemm_cntl_blas,
	                                                 FLA_cntl_flamec_init_i->fla_gemm_cntl_blas,
	                                                 FLA_cntl_flamec_init_i->fla_hemm_cntl_blas,
	                                                 FLA_cntl_flamec_init_i->fla_her2k_cntl_blas );
}

void FLA_Lyap_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flamec_init_i->fla_lyap_cntl_leaf );
	FLA_Cntl_obj_free( FLA_cntl_flamec_init_i->fla_lyap_cntl );

	FLA_Blocksize_free( FLA_cntl_flamec_init_i->fla_lyap_bsize );
}
#endif

extern fla_scal_t*  fla_scal_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_sylv_t*  fla_sylv_cntl;

fla_lyap_t*         fla_lyap_cntl_leaf = NULL;
fla_lyap_t*         fla_lyap_cntl = NULL;
fla_blocksize_t*    fla_lyap_bsize = NULL;

void FLA_Lyap_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_lyap_bsize       = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke an unblocked variant.
	fla_lyap_cntl_leaf   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                                 FLA_UNBLOCKED_VARIANT1,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke a blocked variant.
	fla_lyap_cntl        = FLA_Cntl_lyap_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT1,
	                                                 fla_lyap_bsize,
	                                                 fla_scal_cntl_blas,
	                                                 fla_lyap_cntl_leaf,
	                                                 fla_sylv_cntl,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_hemm_cntl_blas,
	                                                 fla_her2k_cntl_blas );
}

void FLA_Lyap_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_lyap_cntl_leaf );
	FLA_Cntl_obj_free( fla_lyap_cntl );

	FLA_Blocksize_free( fla_lyap_bsize );
}

