/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

fla_swap_t*      fla_swap_cntl_panel = NULL;
fla_swap_t*      fla_swap_cntl_blas = NULL;

fla_tpose_t*     fla_tpose_cntl = NULL;
fla_tpose_t*     fla_tpose_cntl_unb = NULL;
fla_blocksize_t* fla_tpose_bsize = NULL;
fla_blocksize_t* fla_tpose_swap_bsize = NULL;

void FLA_Transpose_cntl_init()
{
	// Set blocksizes based on libgoto query.
	fla_tpose_bsize      = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_tpose_swap_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that performs unblocked variant 2 transposition.
	fla_tpose_cntl_unb   = FLA_Cntl_tpose_obj_create( FLA_FLAT, 
	                                                  FLA_UNBLOCKED_VARIANT2,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that invokes an external implementation of swap.
	fla_swap_cntl_blas   = FLA_Cntl_swap_obj_create( FLA_FLAT,
	                                                 FLA_SUBPROBLEM,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree that invokes unblocked variant 2 of swap.
	fla_swap_cntl_panel  = FLA_Cntl_swap_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT2, 
	                                                 fla_tpose_swap_bsize,
	                                                 fla_swap_cntl_blas );

	// Create a control tree that assumes a large matrix argument.
	fla_tpose_cntl       = FLA_Cntl_tpose_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT2, 
	                                                  fla_tpose_bsize,
	                                                  fla_tpose_cntl_unb,
	                                                  fla_swap_cntl_panel );
}

void FLA_Transpose_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_tpose_cntl );
	FLA_Cntl_obj_free( fla_tpose_cntl_unb );
	FLA_Cntl_obj_free( fla_swap_cntl_panel );
	FLA_Cntl_obj_free( fla_swap_cntl_blas );

	FLA_Blocksize_free( fla_tpose_bsize );
	FLA_Blocksize_free( fla_tpose_swap_bsize );
}

