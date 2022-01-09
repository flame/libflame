/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sident( integer m, float* a, integer a_rs, integer a_cs )
{
	float* alpha;
	integer    i, j;

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

void bl1_dident( integer m, double* a, integer a_rs, integer a_cs )
{
	double* alpha;
	integer     i, j;

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

void bl1_cident( integer m, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* alpha;
	integer       i, j;

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

void bl1_zident( integer m, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* alpha;
	integer       i, j;

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

