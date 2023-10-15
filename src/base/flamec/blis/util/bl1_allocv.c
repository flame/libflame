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

#ifdef BLIS1_ENABLE_USE_OF_FLA_MALLOC
  #include "FLAME.h"
  #define BLIS1_MALLOC FLA_malloc
#else
  #define BLIS1_MALLOC malloc
#endif

void*     bl1_vallocv( uinteger n_elem, uinteger elem_size )
{
	return ( void*  ) BLIS1_MALLOC( n_elem * elem_size );
}

integer*      bl1_iallocv( uinteger n_elem )
{
	return ( integer*   ) BLIS1_MALLOC( n_elem * sizeof( integer ) );
}

float*    bl1_sallocv( uinteger n_elem )
{
	return ( float* ) BLIS1_MALLOC( n_elem * sizeof( float ) );
}

double*   bl1_dallocv( uinteger n_elem )
{
	return ( double* ) BLIS1_MALLOC( n_elem * sizeof( double ) );
}

scomplex* bl1_callocv( uinteger n_elem )
{
	return ( scomplex* ) BLIS1_MALLOC( n_elem * sizeof( scomplex ) );
}

dcomplex* bl1_zallocv( uinteger n_elem )
{
	return ( dcomplex* ) BLIS1_MALLOC( n_elem * sizeof( dcomplex ) );
}

