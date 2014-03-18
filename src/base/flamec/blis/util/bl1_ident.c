
#include "blis1.h"

void bl1_sident( int m, float* a, int a_rs, int a_cs )
{
	float* alpha;
	int    i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0F;

			if ( i == j )
				*alpha = 1.0F;
		}
	}
}

void bl1_dident( int m, double* a, int a_rs, int a_cs )
{
	double* alpha;
	int     i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0;

			if ( i == j )
				*alpha = 1.0;
		}
	}
}

void bl1_cident( int m, scomplex* a, int a_rs, int a_cs )
{
	scomplex* alpha;
	int       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0F;
			alpha->imag = 0.0F;

			if ( i == j )
				alpha->real = 1.0F;
		}
	}
}

void bl1_zident( int m, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* alpha;
	int       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0;
			alpha->imag = 0.0;

			if ( i == j )
				alpha->real = 1.0;
		}
	}
}

