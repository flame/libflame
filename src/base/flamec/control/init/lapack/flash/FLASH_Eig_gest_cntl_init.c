/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_axpy_t*  flash_axpy_cntl;
extern fla_gemm_t*  flash_gemm_cntl_mm;
extern fla_hemm_t*  flash_hemm_cntl_mm;
extern fla_her2k_t* flash_her2k_cntl_mm;
extern fla_trmm_t*  flash_trmm_cntl_mm;
extern fla_trsm_t*  flash_trsm_cntl_mm;

fla_eig_gest_t*     flash_eig_gest_cntl_leaf = NULL;
fla_eig_gest_t*     flash_eig_gest_cntl = NULL;
fla_blocksize_t*    flash_eig_gest_bsize = NULL;

void FLASH_Eig_gest_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_eig_gest_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_eig_gest_cntl_leaf   = FLA_Cntl_eig_gest_obj_create( FLA_HIER,
	                                                           FLA_SUBPROBLEM,
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

	// Create a control tree that assumes A is large.
	flash_eig_gest_cntl        = FLA_Cntl_eig_gest_obj_create( FLA_HIER,
	                                                           FLA_BLOCKED_VARIANT1, 
	                                                           flash_eig_gest_bsize,
	                                                           flash_eig_gest_cntl_leaf,
	                                                           flash_axpy_cntl,
	                                                           flash_axpy_cntl,
	                                                           NULL,
	                                                           NULL,
	                                                           NULL,
	                                                           flash_hemm_cntl_mm,
	                                                           flash_her2k_cntl_mm,
	                                                           flash_trmm_cntl_mm,
	                                                           flash_trmm_cntl_mm,
	                                                           flash_trsm_cntl_mm,
	                                                           flash_trsm_cntl_mm );
}

void FLASH_Eig_gest_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_eig_gest_cntl_leaf );
	FLA_Cntl_obj_free( flash_eig_gest_cntl );

	FLA_Blocksize_free( flash_eig_gest_bsize );
}

