/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_ssymmize( conj1_t conj, uplo1_t uplo, integer m, float* a, integer a_rs, integer a_cs )
{
	float*    a_src;
	float*    a_dst;
	integer       rs_src, cs_src, inc_src;
	integer       rs_dst, cs_dst, inc_dst;
	integer       n_iter;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bl1_scopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );
	}
}

void bl1_dsymmize( conj1_t conj, uplo1_t uplo, integer m, double* a, integer a_rs, integer a_cs )
{
	double*   a_src;
	double*   a_dst;
	integer       rs_src, cs_src, inc_src;
	integer       rs_dst, cs_dst, inc_dst;
	integer       n_iter;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bl1_dcopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );
	}
}

void bl1_csymmize( conj1_t conj, uplo1_t uplo, integer m, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* a_src;
	scomplex* a_dst;
	scomplex* a_jj;
	integer       rs_src, cs_src, inc_src;
	integer       rs_dst, cs_dst, inc_dst;
	integer       n_iter;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bl1_ccopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );

		if ( bl1_is_conj( conj ) )
		{
			a_jj = a + j*a_rs + j*a_cs;
			a_jj->imag = bl1_s0();
		}
	}
}

void bl1_zsymmize( conj1_t conj, uplo1_t uplo, integer m, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* a_src;
	dcomplex* a_dst;
	dcomplex* a_jj;
	integer       rs_src, cs_src, inc_src;
	integer       rs_dst, cs_dst, inc_dst;
	integer       n_iter;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bl1_is_col_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bl1_is_row_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bl1_is_gen_storage( a_rs, a_cs ) && bl1_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bl1_zcopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );

		if ( bl1_is_conj( conj ) )
		{
			a_jj = a + j*a_rs + j*a_cs;
			a_jj->imag = bl1_d0();
		}
	}
}

