
#include "FLAME.h"

FLA_Error FLA_Fill_with_cluster_dist( FLA_Obj n_clusters, FLA_Obj cluster_width, FLA_Obj x )
{
	FLA_Obj      lT,              l0,
	             lB,              l1,
	                              l2;
	FLA_Obj      lT_rest,
	             lT_last;
	FLA_Obj      l, k;
	FLA_Datatype dt_real;
	dim_t        n_x;
	int          nc;
	int          n_regions;
	int          region_width;
	int          leftover_width;
	

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Fill_with_cluster_dist_check( n_clusters, cluster_width, x );

	dt_real = FLA_Obj_datatype_proj_to_real( x );
	n_x     = FLA_Obj_vector_dim( x );

	nc = *FLA_INT_PTR( n_clusters );
	n_regions = 2 * nc;
	region_width   = n_x / n_regions;
	leftover_width = n_x % n_regions;

	// Create a local counter to increment as we create the distribution.
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &k );

	// Create a local vector l. We will work with this vector, which is
	// the same length as x, so that we can use vertical partitioning.
	FLA_Obj_create( dt_real, n_x, 1, 0, 0, &l );

	// Initialize k to 1.
	FLA_Set( FLA_ZERO, k );

	FLA_Part_2x1( l,    &lT,
	                    &lB,            0, FLA_TOP );

	while ( FLA_Obj_length( lT ) < n_regions * region_width )
	{
		FLA_Repart_2x1_to_3x1( lT,                &l0,
		                    /* ** */            /* ******* */
		                                          &l1,
		                       lB,                &l2,        region_width, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		FLA_Fill_with_linear_dist( k, FLA_ONE, l1 );

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &lT,                l0,
		                                              l1,
		                        /* ** */           /* ******* */
		                          &lB,                l2,     FLA_TOP );


		FLA_Part_2x1( lT,   &lT_rest,
		                    &lT_last,            1, FLA_BOTTOM );
		FLA_Copy( lT_last, k );


		FLA_Repart_2x1_to_3x1( lT,                &l0,
		                    /* ** */            /* ******* */
		                                          &l1,
		                       lB,                &l2,        region_width, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		FLA_Fill_with_random_dist( k, cluster_width, l1 );
		FLA_Sort( FLA_FORWARD, l1 );

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &lT,                l0,
		                                              l1,
		                        /* ** */           /* ******* */
		                          &lB,                l2,     FLA_TOP );

		FLA_Part_2x1( lT,   &lT_rest,
		                    &lT_last,            1, FLA_BOTTOM );
		FLA_Copy( lT_last, k );
		FLA_Mult_add( FLA_ONE, FLA_ONE, k );
	}

	if ( leftover_width > 0 )
		FLA_Fill_with_linear_dist( k, FLA_ONE, lB );

	// Normalize by last element.
	//FLA_Part_2x1( l,    &lT,
	//                    &lB,            1, FLA_BOTTOM );
	//FLA_Inv_scal( lB, l );

	// Overwrite x with the distribution we created in l.
	FLA_Copy( l, x );

	FLA_Obj_free( &l );
	FLA_Obj_free( &k );

	return FLA_SUCCESS;
}

