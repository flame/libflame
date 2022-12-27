/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_QR_UT_inc_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* TW )
{
	FLA_Datatype datatype;
	dim_t        m, n;
	dim_t        n_last;
	
	// *** The current QR_UT_inc algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1. We check for that here
	// because we anticipate that we'll use a more general algorithm in the
	// future, and we don't want to forget to remove the constraint. ***
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_QR_UT_inc() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Create hierarchical copy of matrix A_flat.
	FLASH_Obj_create_hier_copy_of_flat( A_flat, depth, b_flash, A );

	// Query the datatype of matrix A_flat.
	datatype = FLA_Obj_datatype( A_flat );
	
	// If the user passed in zero for b_alg, then we need to set the
	// algorithmic (inner) blocksize to a reasonable default value.
	if ( b_alg == 0 )
	{
		b_alg = FLASH_QR_UT_inc_determine_alg_blocksize( *A );
	}

	// Query the element (not scalar) dimensions of the new hierarchical
	// matrix. This is done so we can create T with full blocks for the
	// bottom and right "edge cases" of A.
	m = FLA_Obj_length( *A );
	n = FLA_Obj_width( *A );

	// Create hierarchical matrices T and W. T is lower triangular where
	// each block is b_alg-by-b_flash and W is strictly upper triangular
	// where each block is b_alg-by-b_flash. So we can create them
	// simultaneously as part of the same hierarchical matrix.
	FLASH_Obj_create_ext( datatype, m * b_alg, n * b_flash[0], 
	                      depth, &b_alg, b_flash, 
	                      TW );

	// If the bottom-right-most block along the diagonal is a partial block,
	// adjust the view of the corresponding T block.
	n_last = FLASH_Obj_scalar_width( *A ) % *b_flash;

	if ( n_last > 0 )
	{
		FLA_Obj  TWTL, TWTR,
		         TWBL, TWBR;
		FLA_Obj  TWL,  TWR;
		FLA_Obj  TWT,  TW0,
		         TWB,  TW1,
		               TW2;
		FLA_Obj* TW1p;

		FLA_Part_2x2( *TW,   &TWTL, &TWTR,
		                     &TWBL, &TWBR,    n-1, n-1, FLA_TL );

		FLA_Part_2x1( TWBR,    &TWT,
		                       &TWB,     0, FLA_TOP );

		while ( FLA_Obj_length( TWB ) > 0 )
		{
			FLA_Repart_2x1_to_3x1( TWT,                &TW0, 
			                    /* *** */            /* *** */
			                                           &TW1, 
			                       TWB,                &TW2,        1, FLA_BOTTOM );

			// -----------------------------------------------------------

			TW1p = FLASH_OBJ_PTR_AT( TW1 );

			FLA_Part_1x2( *TW1p,   &TWL, &TWR,    n_last, FLA_LEFT );

			*TW1p = TWL;
			TW1p->m_inner = TW1p->m;
			TW1p->n_inner = TW1p->n;

			// -----------------------------------------------------------

			FLA_Cont_with_3x1_to_2x1( &TWT,                TW0, 
			                                               TW1, 
			                        /* *** */           /* *** */
			                          &TWB,                TW2,     FLA_TOP );
		}


	}
	   
	return FLA_SUCCESS;
}


dim_t FLASH_QR_UT_inc_determine_alg_blocksize( FLA_Obj A )
{
	dim_t b_alg;
	dim_t b_flash;

	// Acquire the storage blocksize.
	b_flash = FLA_Obj_length( *FLASH_OBJ_PTR_AT( A ) );

	// Scale the storage blocksize by a pre-defined scalar to arrive at a
	// reasonable algorithmic blocksize, but make sure it's at least 1.
	b_alg = ( dim_t ) fla_max( ( double ) b_flash * FLA_QR_INNER_TO_OUTER_B_RATIO, 1 );

	return b_alg;
}

