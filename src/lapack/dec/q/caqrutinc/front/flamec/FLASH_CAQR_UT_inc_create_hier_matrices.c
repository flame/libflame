/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_CAQR_UT_inc_create_hier_matrices( dim_t p, FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* ATW, FLA_Obj* R, FLA_Obj* RTW )
{
	FLA_Datatype datatype;
	dim_t        m, n;
	dim_t        nb_part;
	
	// *** The current CAQR_UT_inc algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1.
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_CAQR_UT_inc() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Create hierarchical copy of matrix A_flat.
	FLASH_Obj_create_hier_copy_of_flat( A_flat, depth, b_flash, A );

	// Create hierarchical copy of matrix A_flat.
	FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, *A, R );

	// Query the datatype of matrix A_flat.
	datatype = FLA_Obj_datatype( A_flat );
	
	// If the user passed in zero for b_alg, then we need to set the
	// algorithmic (inner) blocksize to a reasonable default value.
	if ( b_alg == 0 )
	{
		b_alg = FLASH_CAQR_UT_inc_determine_alg_blocksize( *A );
	}

	// Query the element (not scalar) dimensions of the new hierarchical
	// matrix. This is done so we can create T with full blocks for the
	// bottom and right "edge cases" of A.
	m = FLA_Obj_length( *A );
	n = FLA_Obj_width( *A );

	// Create hierarchical matrices T and W for both A and R. T is lower
	// triangular where each block is b_alg-by-b_flash and W is strictly
	// upper triangular where each block is b_alg-by-b_flash. So we can
	// create them simultaneously as part of the same hierarchical matrix.
	FLASH_Obj_create_ext( datatype, m * b_alg, n * b_flash[0], 
	                      depth, &b_alg, b_flash, 
	                      ATW );
	FLASH_Obj_create_ext( datatype, m * b_alg, n * b_flash[0], 
	                      depth, &b_alg, b_flash, 
	                      RTW );

	// If the bottom-right-most block along the diagonal is a partial block,
	// adjust the view of the corresponding T block.
	FLASH_CAQR_UT_inc_adjust_views( *A, *ATW );
	FLASH_CAQR_UT_inc_adjust_views( *A, *RTW );

	// Compute the partition length from the number of partitions.
	nb_part = FLA_CAQR_UT_inc_compute_blocks_per_part( p, *A );

	// Encode block structure (upper tri, full, or zero) into blocks of R.
	FLA_CAQR_UT_inc_init_structure( p, nb_part, *R );

	return FLA_SUCCESS;
}

FLA_Error FLASH_CAQR_UT_inc_adjust_views( FLA_Obj A, FLA_Obj TW )
{
	dim_t b_flash;
	dim_t n, n_last;

	// We can query b_flash as the width of the top-left element of TW.
	b_flash = FLASH_Obj_scalar_width_tl( TW );

	// Query the element (not scalar) n dimension of A.
	n = FLA_Obj_width( A );

	// If the bottom-right-most block along the diagonal is a partial block,
	// adjust the view of the corresponding T block.
	n_last = FLASH_Obj_scalar_width( A ) % b_flash;

	if ( n_last > 0 )
	{
		FLA_Obj  TWTL, TWTR,
		         TWBL, TWBR;
		FLA_Obj  TWL,  TWR;
		FLA_Obj  TWT,  TW0,
		         TWB,  TW1,
		               TW2;
		FLA_Obj* TW1p;

		FLA_Part_2x2( TW,    &TWTL, &TWTR,
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


dim_t FLASH_CAQR_UT_inc_determine_alg_blocksize( FLA_Obj A )
{
	dim_t b_alg;
	dim_t b_flash;

	// Acquire the storage blocksize.
	b_flash = FLA_Obj_length( *FLASH_OBJ_PTR_AT( A ) );

	// Scale the storage blocksize by a pre-defined scalar to arrive at a
	// reasonable algorithmic blocksize, but make sure it's at least 1.
	b_alg = ( dim_t ) max( ( double ) b_flash * FLA_CAQR_INNER_TO_OUTER_B_RATIO, 1 );

	return b_alg;
}

