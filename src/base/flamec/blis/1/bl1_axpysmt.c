/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_saxpysmt( trans1_t trans, integer m, integer n, float* alpha0, float* alpha1, float* a, integer a_rs, integer a_cs, float* beta, float* b, integer b_rs, integer b_cs )
{
	float*    a_begin;
	float*    b_begin;
	float     alpha_prod;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_sscal( n_elem,
		           beta,
		           b_begin, incb );

		bl1_saxpy( n_elem,
		           &alpha_prod,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bl1_daxpysmt( trans1_t trans, integer m, integer n, double* alpha0, double* alpha1, double* a, integer a_rs, integer a_cs, double* beta, double* b, integer b_rs, integer b_cs )
{
	double*   a_begin;
	double*   b_begin;
	double    alpha_prod;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_dscal( n_elem,
		           beta,
		           b_begin, incb );

		bl1_daxpy( n_elem,
		           &alpha_prod,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bl1_caxpysmt( trans1_t trans, integer m, integer n, scomplex* alpha0, scomplex* alpha1, scomplex* a, integer a_rs, integer a_cs, scomplex* beta, scomplex* b, integer b_rs, integer b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	scomplex* a_temp;
	scomplex  alpha_prod;
	integer       inca_temp;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	if ( bl1_does_conj( trans ) )
	{
		conj1_t conj = bl1_proj_trans1_to_conj( trans );

		a_temp = bl1_callocv( n_elem );
		inca_temp = 1;

		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            a_temp,  inca_temp );

			bl1_cscal( n_elem,
			           beta,
			           b_begin, incb );

			bl1_caxpy( n_elem,
			           &alpha_prod,
			           a_temp,  inca_temp, 
			           b_begin, incb );
		}
	
		bl1_cfree( a_temp );
	}
	else // if ( !bl1_does_conj( trans ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_cscal( n_elem,
			           beta,
			           b_begin, incb );

			bl1_caxpy( n_elem,
			           &alpha_prod,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bl1_zaxpysmt( trans1_t trans, integer m, integer n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, integer a_rs, integer a_cs, dcomplex* beta, dcomplex* b, integer b_rs, integer b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	dcomplex* a_temp;
	dcomplex  alpha_prod;
	integer       inca_temp;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	if ( bl1_does_conj( trans ) )
	{
		conj1_t conj = bl1_proj_trans1_to_conj( trans );

		a_temp = bl1_zallocv( n_elem );
		inca_temp = 1;

		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            a_temp,  inca_temp );

			bl1_zscal( n_elem,
			           beta,
			           b_begin, incb );

			bl1_zaxpy( n_elem,
			           &alpha_prod,
			           a_temp,  inca_temp, 
			           b_begin, incb );
		}
	
		bl1_zfree( a_temp );
	}
	else // if ( !bl1_does_conj( trans ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bl1_zscal( n_elem,
			           beta,
			           b_begin, incb );

			bl1_zaxpy( n_elem,
			           &alpha_prod,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

