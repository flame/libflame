/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sshiftdiag( conj1_t conj, integer offset, integer m, integer n, float* sigma, float* a, integer a_rs, integer a_cs )
{
	float* alpha;
	integer    i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha += *sigma;

		++i;
		++j;
	}
}

void bl1_dshiftdiag( conj1_t conj, integer offset, integer m, integer n, double* sigma, double* a, integer a_rs, integer a_cs )
{
	double* alpha;
	integer     i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha += *sigma;

		++i;
		++j;
	}
}

void bl1_csshiftdiag( conj1_t conj, integer offset, integer m, integer n, float* sigma, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* alpha;
	integer       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real += *sigma;

		++i;
		++j;
	}
}

void bl1_zdshiftdiag( conj1_t conj, integer offset, integer m, integer n, double* sigma, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* alpha;
	integer       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real += *sigma;

		++i;
		++j;
	}
}

void bl1_cshiftdiag( conj1_t conj, integer offset, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* alpha;
	scomplex  sigma_conj;
	integer       i, j;

	bl1_ccopys( conj, sigma, &sigma_conj );

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real += sigma_conj.real;
		alpha->imag += sigma_conj.imag;

		++i;
		++j;
	}
}

void bl1_zshiftdiag( conj1_t conj, integer offset, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* alpha;
	dcomplex  sigma_conj;
	integer       i, j;

	bl1_zcopys( conj, sigma, &sigma_conj );

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real += sigma_conj.real;
		alpha->imag += sigma_conj.imag;

		++i;
		++j;
	}
}

