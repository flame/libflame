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

void bl1_abort( void )
{
	abort();
}

void bl1_abort_msg( char* message )
{
	fprintf( stderr, "BLIS: %s\n", message );
	fprintf( stderr, "BLIS: Aborting.\n" );
	bl1_abort();
}

