/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_compute_scaling_ops( integer       m_A,
                                        float*    buff_d, integer inc_d, 
                                        float*    buff_e, integer inc_e,
                                        float*    sigma )
{
	float  one   = bl1_s1();
	float  three = 3.0F;
	float  norm;
	float  eps2;
	float  safmin;
	float  safmax;
	float  ssfmin;
	float  ssfmax;

	// Query some constants.
	eps2   = FLA_Mach_params_ops( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_ops( FLA_MACH_SFMIN );
	safmax = one / safmin;

	// Compute the acceptable range for the 1-norm;
	ssfmax = sqrt( safmax ) / three;
	ssfmin = sqrt( safmin ) / eps2;

	// Compute the 1-norm of the tridiagonal matrix.
	FLA_Norm1_tridiag_ops( m_A,
	                       buff_d, inc_d,
	                       buff_e, inc_e,
	                       &norm );

	// Compute sigma accordingly if norm is outside of the range.
	if ( norm > ssfmax )
	{
		*sigma = ssfmax / norm;
	}
	else if ( norm < ssfmin )
	{
		*sigma = ssfmin / norm;
	}
	else
	{
		*sigma = one;
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_compute_scaling_opd( integer       m_A,
                                        double*   buff_d, integer inc_d, 
                                        double*   buff_e, integer inc_e,
                                        double*   sigma )
{
	double one   = bl1_d1();
	double three = 3.0;
	double norm;
	double eps2;
	double safmin;
	double safmax;
	double ssfmin;
	double ssfmax;

	// Query some constants.
	eps2   = FLA_Mach_params_opd( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );
	safmax = one / safmin;

	// Compute the acceptable range for the 1-norm;
	ssfmax = sqrt( safmax ) / three;
	ssfmin = sqrt( safmin ) / eps2;

	// Compute the 1-norm of the tridiagonal matrix.
	FLA_Norm1_tridiag_opd( m_A,
	                       buff_d, inc_d,
	                       buff_e, inc_e,
	                       &norm );

	// Compute sigma accordingly if norm is outside of the range.
	if ( norm > ssfmax )
	{
		*sigma = ssfmax / norm;
	}
	else if ( norm < ssfmin )
	{
		*sigma = ssfmin / norm;
	}
	else
	{
		*sigma = one;
	}

	return FLA_SUCCESS;
}

