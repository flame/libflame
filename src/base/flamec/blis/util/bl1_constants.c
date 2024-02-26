/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

// --- two ---

float bl1_s2( void )
{
	float x;
	x = 2.0F;
	return x;
}

double bl1_d2( void )
{
	double x;
	x = 2.0;
	return x;
}

scomplex bl1_c2( void )
{
	scomplex x;
	x.real = bl1_s2();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_z2( void )
{
	dcomplex x;
	x.real = bl1_d2();
	x.imag = bl1_d0();
	return x;
}

// --- one ---

float bl1_s1( void )
{
	float x;
	x = 1.0F;
	return x;
}

double bl1_d1( void )
{
	double x;
	x = 1.0;
	return x;
}

scomplex bl1_c1( void )
{
	scomplex x;
	x.real = bl1_s1();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_z1( void )
{
	dcomplex x;
	x.real = bl1_d1();
	x.imag = bl1_d0();
	return x;
}

// --- one half ---

float bl1_s1h( void )
{
	float x;
	x = 0.5F;
	return x;
}

double bl1_d1h( void )
{
	double x;
	x = 0.5;
	return x;
}

scomplex bl1_c1h( void )
{
	scomplex x;
	x.real = bl1_s1h();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_z1h( void )
{
	dcomplex x;
	x.real = bl1_d1h();
	x.imag = bl1_d0();
	return x;
}

// --- zero ---

float bl1_s0( void )
{
	float x;
	x = 0.0F;
	return x;
}

double bl1_d0( void )
{
	double x;
	x = 0.0;
	return x;
}

scomplex bl1_c0( void )
{
	scomplex x;
	x.real = bl1_s0();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_z0( void )
{
	dcomplex x;
	x.real = bl1_d0();
	x.imag = bl1_d0();
	return x;
}

// --- minus one half ---

float bl1_sm1h( void )
{
	float x;
	x = -0.5F;
	return x;
}

double bl1_dm1h( void )
{
	double x;
	x = -0.5;
	return x;
}

scomplex bl1_cm1h( void )
{
	scomplex x;
	x.real = bl1_sm1h();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_zm1h( void )
{
	dcomplex x;
	x.real = bl1_dm1h();
	x.imag = bl1_d0();
	return x;
}

// --- minus one ---

float bl1_sm1( void )
{
	float x;
	x = -1.0F;
	return x;
}

double bl1_dm1( void )
{
	double x;
	x = -1.0;
	return x;
}

scomplex bl1_cm1( void )
{
	scomplex x;
	x.real = bl1_sm1();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_zm1( void )
{
	dcomplex x;
	x.real = bl1_dm1();
	x.imag = bl1_d0();
	return x;
}

// --- minus two ---

float bl1_sm2( void )
{
	float x;
	x = -2.0F;
	return x;
}

double bl1_dm2( void )
{
	double x;
	x = -2.0;
	return x;
}

scomplex bl1_cm2( void )
{
	scomplex x;
	x.real = bl1_sm2();
	x.imag = bl1_s0();
	return x;
}

dcomplex bl1_zm2( void )
{
	dcomplex x;
	x.real = bl1_dm2();
	x.imag = bl1_d0();
	return x;
}

