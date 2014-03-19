/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


FLA_Error FLA_Tevd_find_perfshift_ops( int       m_d,
                                       int       m_l,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e, 
                                       float*    buff_l, int inc_l, 
                                       int*      buff_ls, int inc_ls, 
                                       float*    buff_pu, int inc_pu, 
                                       int*      ij_shift )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_find_perfshift_opd( int       m_d,
                                       int       m_l,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e, 
                                       double*   buff_l, int inc_l, 
                                       int*      buff_ls, int inc_ls, 
                                       double*   buff_pu, int inc_pu, 
                                       int*      ij_shift )
{
	double* d1p;
	double* e1p;
	double* d2p;
	double  wilkshift;
	int     i;
	int     ij_cand;
	double  dist_cand;
	double  pshift_cand;
	
	d1p = buff_d + (m_d-2)*inc_d;
	e1p = buff_e + (m_d-2)*inc_e;
	d2p = buff_d + (m_d-1)*inc_d;

	if ( *buff_ls == -1 )
	{
		*ij_shift = -1;
		return FLA_FAILURE;
	}

	FLA_Wilkshift_tridiag_opd( *d1p,
	                           *e1p,
	                           *d2p,
	                           &wilkshift );

/*
	// If we have shifted here previously, use a Wilkinson shfit.
	prev_shift = buff_pu[ (m_d-1)*inc_pu ];

	if ( prev_shift != 0.0 )
	{
		// *shift = prev_shift;
		*shift = wilkshift;
		return FLA_SUCCESS;
	}
*/
	ij_cand = -1;

	// Find an available (unused) shift.
	for ( i = 0; i < m_l; ++i )
	{
		int* status = buff_ls + (i  )*inc_ls;

		if ( *status == 0 )
		{
			double* lambda1 = buff_l + (i  )*inc_l;
			ij_cand     = i;
			pshift_cand = *lambda1;
			dist_cand   = fabs( wilkshift - pshift_cand );
		}
	}

	if ( ij_cand == -1 )
	{
		*ij_shift = -1;
		*buff_ls  = -1;
		return FLA_FAILURE;
	}

	// Now try to find a shift closer to wilkshift than the
	// first one we found.
	for ( i = 0; i < m_l; ++i )
	{
		double* lambda1 = buff_l  + (i  )*inc_l;
		int*    status  = buff_ls + (i  )*inc_ls;
		double  dist    = fabs( wilkshift - *lambda1 );

		if ( *status == 1 ) continue;

		if ( dist < dist_cand )
		{
			ij_cand = i;
			pshift_cand = *lambda1;
			dist_cand = dist;
		}
	}

	*ij_shift = ij_cand;

	return FLA_SUCCESS;
}

