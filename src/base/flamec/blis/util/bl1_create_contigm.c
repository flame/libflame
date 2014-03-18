
#include "blis1.h"

void bl1_screate_contigm( int m, int n, float* a_save, int a_rs_save, int a_cs_save, float** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bl1_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bl1_sallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bl1_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bl1_scopymt( BLIS1_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bl1_dcreate_contigm( int m, int n, double* a_save, int a_rs_save, int a_cs_save, double** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bl1_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bl1_dallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bl1_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bl1_dcopymt( BLIS1_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bl1_ccreate_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bl1_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bl1_callocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bl1_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bl1_ccopymt( BLIS1_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bl1_zcreate_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bl1_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bl1_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bl1_zallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bl1_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bl1_zcopymt( BLIS1_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

