
#include "blis1.h"

void bl1_sscalm( conj1_t conj, int m, int n, float* alpha, float* a, int a_rs, int a_cs )
{
	float     alpha_conj;
	float*    a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

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

	bl1_scopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_sscal( n_elem,
		           &alpha_conj,
		           a_begin, inca );
	}
}

void bl1_dscalm( conj1_t conj, int m, int n, double* alpha, double* a, int a_rs, int a_cs )
{
	double    alpha_conj;
	double*   a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

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

	bl1_dcopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_dscal( n_elem,
		           &alpha_conj,
		           a_begin, inca );
	}
}

void bl1_csscalm( conj1_t conj, int m, int n, float* alpha, scomplex* a, int a_rs, int a_cs )
{
	float     alpha_conj;
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

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

	bl1_scopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_csscal( n_elem,
		            &alpha_conj,
		            a_begin, inca );
	}
}

void bl1_cscalm( conj1_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs )
{
	scomplex  alpha_conj;
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_ceq1( alpha ) ) return;

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

	bl1_ccopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_cscal( n_elem,
		           &alpha_conj,
		           a_begin, inca );
	}
}

void bl1_zdscalm( conj1_t conj, int m, int n, double* alpha, dcomplex* a, int a_rs, int a_cs )
{
	double    alpha_conj;
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

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

	bl1_dcopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_zdscal( n_elem,
		            &alpha_conj,
		            a_begin, inca );
	}
}

void bl1_zscalm( conj1_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex  alpha_conj;
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_zeq1( alpha ) ) return;

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

	bl1_zcopys( conj, alpha, &alpha_conj );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_zscal( n_elem,
		           &alpha_conj,
		           a_begin, inca );
	}
}

