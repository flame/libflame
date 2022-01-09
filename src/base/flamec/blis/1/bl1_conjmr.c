/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sconjmr( uplo1_t uplo, integer m, integer n, float* a, integer a_rs, integer a_cs )
{
	return;
}

void bl1_dconjmr( uplo1_t uplo, integer m, integer n, double* a, integer a_rs, integer a_cs )
{
	return;
}

void bl1_cconjmr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs )
{
	float   m1 = bl1_sm1();
	float*  a_conj;
	integer     lda, inca;
	integer     n_iter;
	integer     n_elem_max;
	integer     n_elem;
	integer     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		for ( j = 0; j < n_iter; ++j )
		{
			n_elem = bl1_min( j + 1, n_elem_max );
			a_conj = ( float* )( a + j*lda ) + 1;
	
			bl1_sscal( n_elem,
			           &m1,
			           a_conj, 2*inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; ++j )
		{
			n_elem = bl1_max( 0, n_elem_max - j );
			a_conj = ( float* )( a + j*lda + j*inca ) + 1;
	
			if ( n_elem <= 0 ) break;

			bl1_sscal( n_elem,
			           &m1,
			           a_conj, 2*inca );
		}
	}
}

void bl1_zconjmr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs )
{
	double  m1 = bl1_dm1();
	double* a_conj;
	integer     lda, inca;
	integer     n_iter;
	integer     n_elem_max;
	integer     n_elem;
	integer     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		for ( j = 0; j < n_iter; ++j )
		{
			n_elem = bl1_min( j + 1, n_elem_max );
			a_conj = ( double* )( a + j*lda ) + 1;
	
			bl1_dscal( n_elem,
			           &m1,
			           a_conj, 2*inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; ++j )
		{
			n_elem = bl1_max( 0, n_elem_max - j );
			a_conj = ( double* )( a + j*lda + j*inca ) + 1;
	
			if ( n_elem <= 0 ) break;

			bl1_dscal( n_elem,
			           &m1,
			           a_conj, 2*inca );
		}
	}
}
