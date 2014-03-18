
#include "blis1.h"

void bl1_srandm( int m, int n, float* a, int a_rs, int a_cs )
{
	float*    a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_srandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_drandm( int m, int n, double* a, int a_rs, int a_cs )
{
	double*   a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_drandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_crandm( int m, int n, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_crandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_zrandm( int m, int n, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_zrandv( n_elem,
		            a_begin, inca );
	}
}

