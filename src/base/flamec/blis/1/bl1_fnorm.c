/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sfnorm( integer m, integer n, float* a, integer a_rs, integer a_cs, float* norm )
{
	float*    a_ij;
	float     sum;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem;
	integer       i, j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0F;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += (*a_ij) * (*a_ij);
		}
	}
	
	// Compute the norm and store the result.
	*norm = ( float ) sqrt( sum );
}

void bl1_dfnorm( integer m, integer n, double* a, integer a_rs, integer a_cs, double* norm )
{
	double*   a_ij;
	double    sum;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem;
	integer       i, j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += (*a_ij) * (*a_ij);
		}
	}
	
	// Compute the norm and store the result.
	*norm = sqrt( sum );
}

void bl1_cfnorm( integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float* norm )
{
	scomplex* a_ij;
	float     sum;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem;
	integer       i, j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0F;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += a_ij->real * a_ij->real + a_ij->imag * a_ij->imag;
		}
	}
	
	// Compute the norm and store the result.
	*norm = ( float ) sqrt( sum );
}

void bl1_zfnorm( integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double* norm )
{
	dcomplex* a_ij;
	double    sum;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem;
	integer       i, j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += a_ij->real * a_ij->real + a_ij->imag * a_ij->imag;
		}
	}
	
	// Compute the norm and store the result.
	*norm = sqrt( sum );
}

