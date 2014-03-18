
#include "FLAME.h"

FLA_Error FLA_Fill_with_linear_dist( FLA_Obj shift, FLA_Obj delta, FLA_Obj x )
{
	FLA_Obj      lT,              l0,
	             lB,              lambda1,
	                              l2;
	FLA_Obj      l, k, delta2;
	FLA_Datatype dt_real;
	dim_t        n_x;


	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Fill_with_linear_dist_check( shift, delta, x );

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
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &delta2 );
	FLA_Copy( delta, delta2 );

	// Initialize k to shift + delta2.
	FLA_Set( shift, k );
	FLA_Mult_add( FLA_ONE, delta2, k );

	FLA_Part_2x1( l,    &lT,
	                    &lB,            0, FLA_TOP );

	while ( FLA_Obj_length( lB ) > 0 )
	{
		FLA_Repart_2x1_to_3x1( lT,                &l0,
		                    /* ** */            /* ******* */
		                                          &lambda1,
		                       lB,                &l2,        1, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		// lambda1 = k;
		FLA_Copy( k, lambda1 );

		// k = k + delta2;
		FLA_Mult_add( FLA_ONE, delta2, k );

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &lT,                l0,
		                                              lambda1,
		                        /* ** */           /* ******* */
		                          &lB,                l2,     FLA_TOP );
	}

	// Normalize by last element.
	//FLA_Part_2x1( l,    &lT,
	//                    &lB,            1, FLA_BOTTOM );
	//FLA_Inv_scal( lB, l );

	// Overwrite x with the distribution we created in l.
	// If x is complex, then this is where the conversion between
	// datatypes happens.
	FLA_Copy( l, x );

	FLA_Obj_free( &l );
	FLA_Obj_free( &k );
	FLA_Obj_free( &delta2 );

	return FLA_SUCCESS;
}

