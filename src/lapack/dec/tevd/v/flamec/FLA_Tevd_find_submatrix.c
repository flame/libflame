
#include "FLAME.h"


FLA_Error FLA_Tevd_find_submatrix_ops( int       m_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_find_submatrix_opd( int       m_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
	double rzero = bl1_d0();
	double eps;
	int    ij_tl;
	int    ij_br;

	// Initialize some numerical constants.
	eps = FLA_Mach_params_opd( FLA_MACH_EPS );

	// Search for the first non-zero subdiagonal element starting at
	// the index specified by ij_begin.
	for ( ij_tl = ij_begin; ij_tl < m_A - 1; ++ij_tl )
	{
		double* d1     = buff_d + (ij_tl  )*inc_d;
		double* d2     = buff_d + (ij_tl+1)*inc_d;
		double* e1     = buff_e + (ij_tl  )*inc_e;
		double  abs_e1 = fabs( *e1 );

		// If we encounter a non-zero subdiagonal element that is close
		// enough to zero, set it to zero.
		if ( abs_e1 != rzero )
		{
			if ( abs_e1 <= eps * sqrt( fabs( *d1 ) ) *
			                     sqrt( fabs( *d2 ) ) )
			{
#ifdef PRINTF
printf( "FLA_Tevd_find_submatrix_opd: nudging non-zero subdiagonal element (e1) to zero.\n" );
printf( "                             d[%3d]        = %22.19e\n", ij_tl, *d1 );
printf( "                             e[%3d] d[%3d] = %22.19e %22.19e\n", ij_tl, ij_tl+1, *e1, *d2 );
#endif
				*e1 = rzero;
			}
		}

		// If we find a non-zero element, record it and break out of this
		// loop.
		if ( *e1 != rzero )
		{
#ifdef PRINTF
printf( "FLA_Tevd_find_submatrix_opd: found non-zero subdiagonal element\n" );
printf( "                             e[%3d] = %22.19e\n", ij_tl, *e1 );
#endif
			*ijTL = ij_tl;
			break;
		}
	}

	// If ij_tl was incremented all the way up to m_A - 1, then we didn't
	// find any non-zeros.
	if ( ij_tl == m_A - 1 )
	{
#ifdef PRINTF
printf( "FLA_Tevd_find_submatrix_opd: no submatrices found.\n" );
#endif
		return FLA_FAILURE;
	}

	// If we've gotten this far, then a non-zero subdiagonal element was
	// found. Now we must walk the remaining portion of the subdiagonal
	// to find the first zero element, or if one is not found, we simply
	// use the last element of the subdiagonal.
	for ( ij_br = ij_tl; ij_br < m_A - 1; ++ij_br )
	{
		double* d1     = buff_d + (ij_br  )*inc_d;
		double* d2     = buff_d + (ij_br+1)*inc_d;
		double* e1     = buff_e + (ij_br  )*inc_e;
		double  abs_e1 = fabs( *e1 );

		// If we encounter a non-zero subdiagonal element that is close
		// enough to zero, set it to zero.
		if ( abs_e1 != rzero )
		{
			if ( abs_e1 <= eps * sqrt( fabs( *d1 ) ) *
			                     sqrt( fabs( *d2 ) ) )
			{
#ifdef PRINTF
printf( "FLA_Tevd_find_submatrix_opd: nudging non-zero subdiagonal element (e1) to zero.\n" );
printf( "                             d[%3d]        = %22.19e\n", ij_br, *d1 );
printf( "                             e[%3d] d[%3d] = %22.19e %22.19e\n", ij_br, ij_br+1, *e1, *d2 );
#endif
				*e1 = rzero;
			}
		}

		// If we find a zero element, record it and break out of this
		// loop.
		if ( *e1 == rzero )
		{
#ifdef PRINTF
printf( "FLA_Tevd_find_submatrix_opd: found zero subdiagonal element\n" );
printf( "                             e[%3d] = %22.19e\n", ij_br, *e1 );
#endif
			break;
		}
	}

	// If a zero element was found, then ij_br should hold the index of
	// that element. If a zero element was not found, then ij_br should
	// hold m_A - 1. Either way, we save the value and return success.
	*ijBR = ij_br;

	return FLA_SUCCESS;
}

