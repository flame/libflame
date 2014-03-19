/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_smaxabsv( int n, float* x, int incx, float* maxabs )
{
	float*    chi;
	float     maxabs_cand;
	float     maxabs_temp;
	int       i;

	bl1_sabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_sabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bl1_dmaxabsv( int n, double* x, int incx, double* maxabs )
{
	double*   chi;
	double    maxabs_cand;
	double    maxabs_temp;
	int       i;

	bl1_dabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_dabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bl1_cmaxabsv( int n, scomplex* x, int incx, float* maxabs )
{
	scomplex* chi;
	float     maxabs_cand;
	float     maxabs_temp;
	int       i;

	bl1_csabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_csabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bl1_zmaxabsv( int n, dcomplex* x, int incx, double* maxabs )
{
	dcomplex* chi;
	double    maxabs_cand;
	double    maxabs_temp;
	int       i;

	bl1_zdabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_zdabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}


