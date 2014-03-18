
#include "blis1.h"

void bl1_set_contig_strides( int m, int n, int* rs, int* cs )
{
	// Default to column-major order.
	*rs = 1;
	*cs = m;

	// Handle special cases first.
	// Check the strides, and modify them if needed.
	if ( *rs == 1 && *cs == 1 )
	{
		// If both strides are unit, we are probably trying to create a
		// 1-by-n matrix in column-major order, or an m-by-1 matrix in
		// row-major order. We have decided to "reserve" the case where
		// rs == cs == 1 for scalars only, as having unit strides can
		// upset the BLAS error checking when attempting to induce a
		// row-major operation.
		if ( m > 1 && n == 1 )
		{
			// Set the column stride to indicate that this is an m-by-1
			// matrix (or vector) stored in column-major order. This is
			// necessary because, in some cases, we have to satisfy error
			// checking in the underlying BLAS library, which expects the
			// leading dimension to be set to at least m, even if it will
			// never be used for indexing since there is only one column
			// of data. Note that rs is already set to 1.
			*cs = m;
		}
		else if ( m == 1 && 1 < n )
		{
			// Set the row stride to indicate that this is a 1-by-n matrix
			// stored in row-major order. Note that cs is already set to 1.
			*rs = n;
		}
		else
		{
			// If m == n == 1, then we are dealing with a scalar. Since rs
			// and cs do not exceed m and n, we don't have to do anything.
		}
	}
}

