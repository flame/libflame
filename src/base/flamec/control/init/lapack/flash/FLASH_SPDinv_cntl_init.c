/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLASH_SPDinv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	// Rather than embed a blocksize, we store the cutoff matrix size for
	// switching from external routines to internal FLAME variants.
	FLA_cntl_flash_init_i->flash_spdinv_size_cutoff = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Initialize a control tree node that calls the top-level Cholesky
	// factorization, Trinagular inversion, and Triangular-transpose matrix
	// multiply control trees. 
	FLA_cntl_flash_init_i->flash_spdinv_cntl        = FLA_Cntl_spdinv_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT1, 
	                                                       FLA_cntl_flash_init_i->flash_spdinv_size_cutoff,
	                                                       FLA_cntl_flash_init_i->flash_chol_cntl,
	                                                       FLA_cntl_flash_init_i->flash_trinv_cntl,
	                                                       FLA_cntl_flash_init_i->flash_ttmm_cntl );
}

void FLASH_SPDinv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_spdinv_cntl );

	FLA_Blocksize_free( FLA_cntl_flash_init_i->flash_spdinv_size_cutoff );
}

#endif

extern fla_chol_t*  flash_chol_cntl;
extern fla_trinv_t* flash_trinv_cntl;
extern fla_ttmm_t*  flash_ttmm_cntl;

fla_spdinv_t*       flash_spdinv_cntl = NULL;
fla_blocksize_t*    flash_spdinv_size_cutoff = NULL;

void FLASH_SPDinv_cntl_init()
{
	// Rather than embed a blocksize, we store the cutoff matrix size for
	// switching from external routines to internal FLAME variants.
	flash_spdinv_size_cutoff = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Initialize a control tree node that calls the top-level Cholesky
	// factorization, Trinagular inversion, and Triangular-transpose matrix
	// multiply control trees. 
	flash_spdinv_cntl        = FLA_Cntl_spdinv_obj_create( FLA_HIER,
	                                                       FLA_BLOCKED_VARIANT1, 
	                                                       flash_spdinv_size_cutoff,
	                                                       flash_chol_cntl,
	                                                       flash_trinv_cntl,
	                                                       flash_ttmm_cntl );
}

void FLASH_SPDinv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_spdinv_cntl );

	FLA_Blocksize_free( flash_spdinv_size_cutoff );
}

