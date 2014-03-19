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

void*     bl1_vallocv( unsigned int n_elem, unsigned int elem_size )
{
	return ( void*  ) BLIS1_MALLOC( n_elem * elem_size );
}

int*      bl1_iallocv( unsigned int n_elem )
{
	return ( int*   ) BLIS1_MALLOC( n_elem * sizeof( int ) );
}

float*    bl1_sallocv( unsigned int n_elem )
{
	return ( float* ) BLIS1_MALLOC( n_elem * sizeof( float ) );
}

double*   bl1_dallocv( unsigned int n_elem )
{
	return ( double* ) BLIS1_MALLOC( n_elem * sizeof( double ) );
}

scomplex* bl1_callocv( unsigned int n_elem )
{
	return ( scomplex* ) BLIS1_MALLOC( n_elem * sizeof( scomplex ) );
}

dcomplex* bl1_zallocv( unsigned int n_elem )
{
	return ( dcomplex* ) BLIS1_MALLOC( n_elem * sizeof( dcomplex ) );
}

