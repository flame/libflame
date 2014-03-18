
#include "blis1.h"

void bl1_sscalmr( uplo1_t uplo, int m, int n, float* alpha, float* a, int a_rs, int a_cs )
{
	float*    a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_sscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_sscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
}

void bl1_dscalmr( uplo1_t uplo, int m, int n, double* alpha, double* a, int a_rs, int a_cs )
{
	double*   a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_dscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_dscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
}

void bl1_csscalmr( uplo1_t uplo, int m, int n, float* alpha, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_csscal( n_elem,
			            alpha,
			            a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_csscal( n_elem,
			            alpha,
			            a_begin, inca );
		}
	}
}

void bl1_cscalmr( uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_ceq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_cscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_cscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
}

void bl1_zdscalmr( uplo1_t uplo, int m, int n, double* alpha, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_zdscal( n_elem,
			            alpha,
			            a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_zdscal( n_elem,
			            alpha,
			            a_begin, inca );
		}
	}
}

void bl1_zscalmr( uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;
	if ( bl1_zeq1( alpha ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix
	// by rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_zscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			if ( n_elem <= 0 ) break;

			bl1_zscal( n_elem,
			           alpha,
			           a_begin, inca );
		}
	}
}

