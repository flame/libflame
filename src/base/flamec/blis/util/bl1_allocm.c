/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

#ifdef BLIS1_ENABLE_USE_OF_FLA_MALLOC
  #include "FLAME.h"
  #define BLIS1_MALLOC FLA_malloc
#else
  #define BLIS1_MALLOC malloc
#endif

void*     bl1_vallocm( unsigned int m, unsigned int n, unsigned int elem_size )
{
	return ( void* ) BLIS1_MALLOC( m * n * elem_size );
}

int*      bl1_iallocm( unsigned int m, unsigned int n )
{
	return ( int* ) BLIS1_MALLOC( m * n * sizeof( int ) );
}

float*    bl1_sallocm( unsigned int m, unsigned int n )
{
	return ( float* ) BLIS1_MALLOC( m * n * sizeof( float ) );
}

double*   bl1_dallocm( unsigned int m, unsigned int n )
{
	return ( double* ) BLIS1_MALLOC( m * n * sizeof( double ) );
}

scomplex* bl1_callocm( unsigned int m, unsigned int n )
{
	return ( scomplex* ) BLIS1_MALLOC( m * n * sizeof( scomplex ) );
}

dcomplex* bl1_zallocm( unsigned int m, unsigned int n )
{
	return ( dcomplex* ) BLIS1_MALLOC( m * n * sizeof( dcomplex ) );
}

