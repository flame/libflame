/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_srandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, float* a, integer a_rs, integer a_cs )
{
	float*    a_begin;
	float*    ajj;
	float     one;
	float     zero;
	float     ord;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize some scalars.
	one      = bl1_s1();
	zero     = bl1_s0();
	ord      = ( float ) bl1_max( m, n );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bl1_srandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bl1_sinvscalv( BLIS1_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_srands( ajj );
					bl1_sabsval2( ajj, ajj );
					bl1_sadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bl1_ssetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bl1_ssetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_srands( ajj );
					bl1_sabsval2( ajj, ajj );
					bl1_sadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bl1_srandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bl1_sinvscalv( BLIS1_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bl1_drandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, double* a, integer a_rs, integer a_cs )
{
	double*   a_begin;
	double*   ajj;
	double    one;
	double    zero;
	double    ord;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize some scalars.
	one      = bl1_d1();
	zero     = bl1_d0();
	ord      = ( double ) bl1_max( m, n );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bl1_drandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bl1_dinvscalv( BLIS1_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_drands( ajj );
					bl1_dabsval2( ajj, ajj );
					bl1_dadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bl1_dsetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bl1_dsetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_drands( ajj );
					bl1_dabsval2( ajj, ajj );
					bl1_dadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bl1_drandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bl1_dinvscalv( BLIS1_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bl1_crandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* a_begin;
	scomplex* ajj;
	scomplex  one;
	scomplex  zero;
	scomplex  ord;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize some scalars.
	one      = bl1_c1();
	zero     = bl1_c0();
	ord      = bl1_c0();
	ord.real = ( float ) bl1_max( m, n );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bl1_crandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bl1_cinvscalv( BLIS1_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_crands( ajj );
					bl1_cabsval2( ajj, ajj );
					bl1_cadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bl1_csetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bl1_csetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_crands( ajj );
					bl1_cabsval2( ajj, ajj );
					bl1_cadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bl1_crandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bl1_cinvscalv( BLIS1_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bl1_zrandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* a_begin;
	dcomplex* ajj;
	dcomplex  one;
	dcomplex  zero;
	dcomplex  ord;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize some scalars.
	one      = bl1_z1();
	zero     = bl1_z0();
	ord      = bl1_z0();
	ord.real = ( double ) bl1_max( m, n );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bl1_zrandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bl1_zinvscalv( BLIS1_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_zrands( ajj );
					bl1_zabsval2( ajj, ajj );
					bl1_zadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bl1_zsetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bl1_zsetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bl1_is_unit_diag( diag ) )    *ajj = one;
				else if ( bl1_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bl1_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bl1_zrands( ajj );
					bl1_zabsval2( ajj, ajj );
					bl1_zadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bl1_zrandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bl1_zinvscalv( BLIS1_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

