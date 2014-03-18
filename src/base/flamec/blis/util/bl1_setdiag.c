
#include "blis1.h"

void bl1_isetdiag( int offset, int m, int n, int* sigma, int* a, int a_rs, int a_cs )
{
	int*   alpha;
	int    i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha = *sigma;

		++i;
		++j;
	}
}

void bl1_ssetdiag( int offset, int m, int n, float* sigma, float* a, int a_rs, int a_cs )
{
	float* alpha;
	int    i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha = *sigma;

		++i;
		++j;
	}
}

void bl1_dsetdiag( int offset, int m, int n, double* sigma, double* a, int a_rs, int a_cs )
{
	double* alpha;
	int     i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha = *sigma;

		++i;
		++j;
	}
}

void bl1_csetdiag( int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs )
{
	scomplex* alpha;
	int       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real = sigma->real;
		alpha->imag = sigma->imag;

		++i;
		++j;
	}
}

void bl1_zsetdiag( int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* alpha;
	int       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real = sigma->real;
		alpha->imag = sigma->imag;

		++i;
		++j;
	}
}

