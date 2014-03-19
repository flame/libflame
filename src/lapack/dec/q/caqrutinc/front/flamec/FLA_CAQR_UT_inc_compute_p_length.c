/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

dim_t FLA_CAQR_UT_inc_compute_blocks_per_part( dim_t p, FLA_Obj A )
{
	dim_t nb_part;
	dim_t nb_left;
	dim_t num_blocks;

	// Query the element (not scalar) length of A.
	num_blocks = FLA_Obj_length( A );
	
	// Compute the number of blocks per partitions.
	nb_part = num_blocks / p;
	nb_left = num_blocks % p;

	// If there are leftover blocks, increase nb_part by one.
	if ( nb_left > 0 ) nb_part += 1;

	return nb_part;
}

