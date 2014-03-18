
#include "blis1.h"

void bl1_sshiftdiag( conj1_t conj, int offset, int m, int n, float* sigma, float* a, int a_rs, int a_cs )
{
	float* alpha;
	int    i, j;

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

void bl1_dshiftdiag( conj1_t conj, int offset, int m, int n, double* sigma, double* a, int a_rs, int a_cs )
{
	double* alpha;
	int     i, j;

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

void bl1_csshiftdiag( conj1_t conj, int offset, int m, int n, float* sigma, scomplex* a, int a_rs, int a_cs )
{
	scomplex* alpha;
	int       i, j;

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

void bl1_zdshiftdiag( conj1_t conj, int offset, int m, int n, double* sigma, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* alpha;
	int       i, j;

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

void bl1_cshiftdiag( conj1_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs )
{
	scomplex* alpha;
	scomplex  sigma_conj;
	int       i, j;

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

void bl1_zshiftdiag( conj1_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* alpha;
	dcomplex  sigma_conj;
	int       i, j;

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

