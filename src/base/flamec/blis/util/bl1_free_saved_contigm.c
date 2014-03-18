
#include "blis1.h"

void bl1_sfree_saved_contigm( int m, int n, float* a_save, int a_rs_save, int a_cs_save, float** a, int* a_rs, int* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bl1_scopymt( BLIS1_NO_TRANSPOSE,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bl1_sfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_dfree_saved_contigm( int m, int n, double* a_save, int a_rs_save, int a_cs_save, double** a, int* a_rs, int* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bl1_dcopymt( BLIS1_NO_TRANSPOSE,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bl1_dfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_cfree_saved_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bl1_ccopymt( BLIS1_NO_TRANSPOSE,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bl1_cfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bl1_zfree_saved_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs )
{
	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bl1_zcopymt( BLIS1_NO_TRANSPOSE,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bl1_zfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

