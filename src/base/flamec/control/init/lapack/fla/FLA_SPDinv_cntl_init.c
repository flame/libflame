/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_chol_t*  fla_chol_cntl;
extern fla_trinv_t* fla_trinv_cntl;
extern fla_ttmm_t*  fla_ttmm_cntl;

fla_spdinv_t*       fla_spdinv_cntl = NULL;
fla_blocksize_t*    fla_spdinv_size_cutoff = NULL;

void FLA_SPDinv_cntl_init()
{
	// Rather than embed a blocksize, we store the cutoff matrix size for
	// switching from external routines to internal FLAME variants.
	fla_spdinv_size_cutoff = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Initialize a control tree node that calls the top-level Cholesky
	// factorization, Trinagular inversion, and Triangular-transpose matrix
	// multiply control trees. 
	fla_spdinv_cntl        = FLA_Cntl_spdinv_obj_create( FLA_FLAT,
	                                                     FLA_BLOCKED_VARIANT1, 
	                                                     fla_spdinv_size_cutoff,
	                                                     fla_chol_cntl,
	                                                     fla_trinv_cntl,
	                                                     fla_ttmm_cntl );
}

void FLA_SPDinv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_spdinv_cntl );

	FLA_Blocksize_free( fla_spdinv_size_cutoff );
}

