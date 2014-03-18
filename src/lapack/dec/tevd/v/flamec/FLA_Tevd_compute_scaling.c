
#include "FLAME.h"

FLA_Error FLA_Tevd_compute_scaling_ops( int       m_A,
                                        float*    buff_d, int inc_d, 
                                        float*    buff_e, int inc_e,
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

FLA_Error FLA_Tevd_compute_scaling_opd( int       m_A,
                                        double*   buff_d, int inc_d, 
                                        double*   buff_e, int inc_e,
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

