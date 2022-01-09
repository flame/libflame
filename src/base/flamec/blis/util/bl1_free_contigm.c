/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sfree_contigm( float* a_save, integer a_rs_save, integer a_cs_save, float** a, integer* a_rs, integer* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Free the temporary contiguous storage for the matrix.
		bl1_sfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_dfree_contigm( double* a_save, integer a_rs_save, integer a_cs_save, double** a, integer* a_rs, integer* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Free the temporary contiguous storage for the matrix.
		bl1_dfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_cfree_contigm( scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Free the temporary contiguous storage for the matrix.
		bl1_cfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_zfree_contigm( dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Free the temporary contiguous storage for the matrix.
		bl1_zfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

