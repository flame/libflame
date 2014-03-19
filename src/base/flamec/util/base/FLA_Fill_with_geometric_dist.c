/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fill_with_geometric_dist( FLA_Obj alpha, FLA_Obj x )
{
	FLA_Obj      lT,              l0,
	             lB,              lambda1,
	                              l2;
	FLA_Obj      l, k, alpha2, temp;
	FLA_Datatype dt_real;
	dim_t        n_x;


	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Fill_with_geometric_dist_check( alpha, x );

	dt_real = FLA_Obj_datatype_proj_to_real( x );
	n_x     = FLA_Obj_vector_dim( x );

	// Create a local counter to increment as we create the distribution.
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &k );

	// Create a local vector l. We will work with this vector, which is
	// the same length as x, so that we can use vertical partitioning.
	FLA_Obj_create( dt_real, n_x, 1, 0, 0, &l );

	// Create a local real scalar alpha2 of the same precision as
	// alpha. Then copy alpha to alpha2, which will convert the
	// complex value to real, if necessary (ie: if alpha is complex).
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &alpha2 );
	FLA_Copy( alpha, alpha2 );

	// Create a temporary scalar.
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &temp );

	// Initialize k to 0.
	FLA_Set( FLA_ZERO, k );

	FLA_Part_2x1( l,    &lT,
	                    &lB,            0, FLA_TOP );

	while ( FLA_Obj_length( lB ) > 0 )
	{
		FLA_Repart_2x1_to_3x1( lT,                &l0,
		                    /* ** */            /* ******* */
		                                          &lambda1,
		                       lB,                &l2,        1, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		// lambda1 = alpha * (1 - alpha)^k;
		FLA_Set( FLA_ONE, temp );
		FLA_Mult_add( FLA_MINUS_ONE, alpha2, temp );
		FLA_Pow( temp, k, lambda1 );
		FLA_Scal( alpha2, lambda1 );

		// k = k + 1;
		FLA_Mult_add( FLA_ONE, FLA_ONE, k );

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &lT,                l0,
		                                              lambda1,
		                        /* ** */           /* ******* */
		                          &lB,                l2,     FLA_TOP );
	}

	// Normalize by first element.
	//FLA_Part_2x1( l,    &lT,
	//                    &lB,            1, FLA_TOP );
	//FLA_Inv_scal( lT, l );

	// Overwrite x with the distribution we created in l.
	// If x is complex, then this is where the conversion between
	// datatypes happens.
	FLA_Copy( l, x );

	FLA_Obj_free( &l );
	FLA_Obj_free( &k );
	FLA_Obj_free( &alpha2 );
	FLA_Obj_free( &temp );

	return FLA_SUCCESS;
}

