
#include "blis1.h"

void bl1_scopymr( uplo1_t uplo, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) && bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_scopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_scopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bl1_dcopymr( uplo1_t uplo, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) && bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_dcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_dcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bl1_ccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) && bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_ccopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_ccopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bl1_zcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) && bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_zcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

// --- Mixed-datatype and general stride copy routines---------------

// ss
void bl1_sscopymr( uplo1_t uplo, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_scopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_scopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// sd ds
void bl1_sdcopymr( uplo1_t uplo, int m, int n, float* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	float*    a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_sdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_sdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_dscopymr( uplo1_t uplo, int m, int n, double* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	double*   a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_dscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_dscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// sc cs
void bl1_sccopymr( uplo1_t uplo, int m, int n, float* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_sccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_sccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_cscopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_cscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_cscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// sz zs
void bl1_szcopymr( uplo1_t uplo, int m, int n, float* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_szcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_szcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_zscopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_zscopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// dd
void bl1_ddcopymr( uplo1_t uplo, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_dcopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_dcopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// dc cd
void bl1_dccopymr( uplo1_t uplo, int m, int n, double* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_dccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_dccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_cdcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_cdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_cdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// dz zd
void bl1_dzcopymr( uplo1_t uplo, int m, int n, double* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_dzcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_dzcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_zdcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_zdcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// cc
void bl1_cccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_ccopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_ccopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// cz zc
void bl1_czcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_czcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_czcopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bl1_zccopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_zccopyv( BLIS1_NO_CONJUGATE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// zz
void bl1_zzcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_swap_ints( ldb, incb );
		bl1_toggle_uplo( uplo );
	}
	
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zcopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bl1_zcopyv( BLIS1_NO_CONJUGATE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

