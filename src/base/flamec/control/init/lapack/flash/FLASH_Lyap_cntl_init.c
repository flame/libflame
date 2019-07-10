/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_scal_t*  flash_scal_cntl;
extern fla_gemm_t*  flash_gemm_cntl_pm;
extern fla_hemm_t*  flash_hemm_cntl_mp;
extern fla_her2k_t* flash_her2k_cntl_ip;

extern fla_sylv_t*  flash_sylv_cntl;

fla_lyap_t*         flash_lyap_cntl_leaf = NULL;
fla_lyap_t*         flash_lyap_cntl = NULL;
fla_blocksize_t*    flash_lyap_bsize = NULL;

void FLASH_Lyap_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_lyap_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and C are b x b blocks.
	flash_lyap_cntl_leaf   = FLA_Cntl_lyap_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is a matrix and C is a matrix.
	flash_lyap_cntl        = FLA_Cntl_lyap_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   flash_lyap_bsize,
	                                                   flash_scal_cntl,
	                                                   flash_lyap_cntl_leaf,
	                                                   flash_sylv_cntl,
	                                                   NULL, //flash_gemm_cntl_pm,
	                                                   NULL, //flash_gemm_cntl_pm,
	                                                   flash_hemm_cntl_mp,
	                                                   flash_her2k_cntl_ip );
}

void FLASH_Lyap_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_lyap_cntl_leaf );
	FLA_Cntl_obj_free( flash_lyap_cntl );

	FLA_Blocksize_free( flash_lyap_bsize );
}

