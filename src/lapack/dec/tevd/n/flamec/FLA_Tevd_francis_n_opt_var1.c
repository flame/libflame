/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_francis_n_opt_var1( FLA_Obj shift, FLA_Obj d, FLA_Obj e )
{
	FLA_Datatype datatype;
	integer          m_A;
	integer          inc_d;
	integer          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_shift = FLA_FLOAT_PTR( shift );
			float*    buff_d     = FLA_FLOAT_PTR( d );
			float*    buff_e     = FLA_FLOAT_PTR( e );

			FLA_Tevd_francis_n_ops_var1( m_A,
			                             buff_shift,
			                             buff_d, inc_d,
			                             buff_e, inc_e );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_shift = FLA_DOUBLE_PTR( shift );
			double*   buff_d     = FLA_DOUBLE_PTR( d );
			double*   buff_e     = FLA_DOUBLE_PTR( e );

			FLA_Tevd_francis_n_opd_var1( m_A,
			                             buff_shift,
			                             buff_d, inc_d,
			                             buff_e, inc_e );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_francis_n_ops_var1( integer       m_A,
                                       float*    buff_shift,
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e )
{
	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_francis_n_opd_var1( integer       m_A,
                                       double*   buff_shift,
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e )
{
	double    eps2, safmin;
	double    temp0, temp1;
	double    bulge;
	double    gamma, sigma;
	integer       ij_deflated;
	integer       i;

	// Initialize the deflation index.
	ij_deflated = FLA_SUCCESS;

	// Initialize the bulge variable to zero.
	bulge = 0.0;

	// Query epsilon and safmin.
	eps2   = FLA_Mach_params_opd( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );

	// Apply the rotations in forward order.
	for ( i = 0; i < m_A - 1; ++i )
	{
		double*   alpha00  = buff_d + (i-1)*inc_d;
		double*   alpha10  = buff_e + (i-1)*inc_e;
		double*   alpha20  = &bulge;

		double*   alpha11  = buff_d + (i  )*inc_d;
		double*   alpha21  = buff_e + (i  )*inc_e;
		double*   alpha22  = buff_d + (i+1)*inc_d;

		double*   alpha31  = &bulge;
		double*   alpha32  = buff_e + (i+1)*inc_e;

		double*   gamma1   = &gamma;
		double*   sigma1   = &sigma;

		double    alpha10_new;

		integer       m_behind = i;
		integer       m_ahead  = m_A - i - 2;

		/*------------------------------------------------------------*/

		if ( i == 0 )
		{
			// Induce an iteration that introduces the bulge by
			// changing the addresses of alpha10 and alpha20.
			temp0 = *buff_d - *buff_shift;
			temp1 = *buff_e;
			alpha10 = &temp0;
			alpha20 = &temp1;

			// Compute a new Givens rotation that introduces the bulge.
			MAC_Givens2_opd( alpha10,
			                 alpha20,
			                 gamma1,
			                 sigma1,
			                 &alpha10_new );

			// We don't apply the Givens rotation to the 2x1 vector at
			// alpha10 when introducing the bulge.
		}
		else
		{
			// Compute a new Givens rotation to push the bulge.
			MAC_Givens2_opd( alpha10,
			                 alpha20,
			                 gamma1,
			                 sigma1,
			                 &alpha10_new );

			// Apply the Givens rotation to the 2x1 vector from which it
			// was computed, which annihilates alpha20.
			*alpha10 = alpha10_new;
			*alpha20 = 0.0;
		}

		// Apply the Givens rotation to the 2x2 submatrix at alpha11.
		MAC_Apply_GTG_opd( gamma1,
		                   sigma1,
		                   alpha11,
		                   alpha21,
		                   alpha22 );

		if ( m_ahead > 0 )
		{
			// Apply the Givens rotation to the 1x2 vector below the 2x2
			// submatrix. This should move the bulge to alpha31.
			MAC_Apply_G_1x2_opd( gamma1,
			                     sigma1,
			                     alpha31,
			                     alpha32 );

			// Check for deflation after applying the rotations, except after
			// applying the initial bulge-introducing rotations.
			if ( m_behind > 0 )
			{
				// We check for deflation in the previous column now that we
				// are done modifying it. If deflation occurred, record the
				// index.
				if ( MAC_Tevd_eigval_converged2_opd( eps2, safmin, *alpha00, *alpha10, *alpha11 ) )
				{
					ij_deflated = i - 1;
				}
			}

			// Sanity check. If the bulge ever disappears, something is wrong.
			if ( *alpha31 == 0.0 )
			{
				printf( "FLA_Tevd_francis_n_opt_var1: bulge disappeared!\n" );
				if ( MAC_Tevd_eigval_converged2_opd( eps2, safmin, *alpha11, *alpha21, *alpha22 ) )
				{
					printf( "FLA_Tevd_francis_n_opt_var1: deflation detected (col %d)\n", (int) i );
					printf( "FLA_Tevd_francis_n_opt_var1: alpha11         = %23.19e\n", *alpha11 );
					printf( "FLA_Tevd_francis_n_opt_var1: alpha21 alpha22 = %23.19e %23.19e\n", *alpha21, *alpha22 );
					return i;
				}
				else
				{
					printf( "FLA_Tevd_francis_n_opt_var1: but NO deflation detected! (col %d)\n", (int) i );
					printf( "FLA_Tevd_francis_n_opt_var1: alpha11         = %23.19e\n", *alpha11 );
					printf( "FLA_Tevd_francis_n_opt_var1: alpha21 alpha22 = %23.19e %23.19e\n", *alpha21, *alpha22 );
					FLA_Abort();
					return FLA_FAILURE;
				}
			}
		}

		/*------------------------------------------------------------*/
	}

	// Return the index of column where deflation most recently occurred,
	// or FLA_SUCCESS if no deflation was detected.
	return ij_deflated;
}

