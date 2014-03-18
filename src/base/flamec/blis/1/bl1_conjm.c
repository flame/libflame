
#include "blis1.h"

void bl1_sconjm( int m, int n, float* a, int a_rs, int a_cs )
{
	return;
}

void bl1_dconjm( int m, int n, double* a, int a_rs, int a_cs )
{
	return;
}

void bl1_cconjm( int m, int n, scomplex* a, int a_rs, int a_cs )
{
	float   m1 = bl1_sm1();
	float*  a_conj;
	int     lda, inca;
	int     n_iter;
	int     n_elem;
	int     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
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

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	for ( j = 0; j < n_iter; ++j )
	{
		a_conj = ( float* )( a + j*lda ) + 1;

		bl1_sscal( n_elem,
		           &m1,
		           a_conj, 2*inca );
	}
}

void bl1_zconjm( int m, int n, dcomplex* a, int a_rs, int a_cs )
{
	double  m1 = bl1_dm1();
	double* a_conj;
	int     lda, inca;
	int     n_iter;
	int     n_elem;
	int     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
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

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	for ( j = 0; j < n_iter; ++j )
	{
		a_conj = ( double* )( a + j*lda ) + 1;

		bl1_dscal( n_elem,
		           &m1,
		           a_conj, 2*inca );
	}
}
