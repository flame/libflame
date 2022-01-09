/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_opt_var2( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of G:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  This variant iterates in pipelined, overlapping fashion and
  applies rotations to two columns at a time.

  -FGVZ
*/
{
	FLA_Datatype datatype;
	integer          k_G, m_A, n_A;
	integer          rs_G, cs_G;
	integer          rs_A, cs_A;

	datatype = FLA_Obj_datatype( A );

	k_G      = FLA_Obj_width( G );
	m_A      = FLA_Obj_length( A );
	n_A      = FLA_Obj_width( A );

	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_A     = FLA_Obj_row_stride( A );
	cs_A     = FLA_Obj_col_stride( A );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			float*    buff_A = ( float*    ) FLA_FLOAT_PTR( A );

			FLA_Apply_G_rf_ops_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

			FLA_Apply_G_rf_opd_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			FLA_Apply_G_rf_opc_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			FLA_Apply_G_rf_opz_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_ops_var2( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   scomplex* buff_G, integer rs_G, integer cs_G,
                                   float*    buff_A, integer rs_A, integer cs_A )
{
	float     one  = bl1_s1();
	float     zero = bl1_s0();
	float     gamma;
	float     sigma;
	float*    a1;
	float*    a2;
	scomplex* g11;
	integer       j, g, k;
	integer       nG, nG_app;
	integer       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_ops_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = 0; j < k_minus_1; ++j )
	{
		nG_app = j + 1;

		for ( k = 0, g = nG_app - 1; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ops( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ops( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ops( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opd_var2( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   dcomplex* buff_G, integer rs_G, integer cs_G,
                                   double*   buff_A, integer rs_A, integer cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma;
	double    sigma;
	double*   a1;
	double*   a2;
	dcomplex* g11;
	integer       j, g, k;
	integer       nG, nG_app;
	integer       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_opd_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = 0; j < k_minus_1; ++j )
	{
		nG_app = j + 1;

		for ( k = 0, g = nG_app - 1; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opc_var2( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   scomplex* buff_G, integer rs_G, integer cs_G,
                                   scomplex* buff_A, integer rs_A, integer cs_A )
{
	float     one  = bl1_s1();
	float     zero = bl1_s0();
	float     gamma;
	float     sigma;
	scomplex* a1;
	scomplex* a2;
	scomplex* g11;
	integer       j, g, k;
	integer       nG, nG_app;
	integer       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_opc_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = 0; j < k_minus_1; ++j )
	{
		nG_app = j + 1;

		for ( k = 0, g = nG_app - 1; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opz_var2( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   dcomplex* buff_G, integer rs_G, integer cs_G,
                                   dcomplex* buff_A, integer rs_A, integer cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma;
	double    sigma;
	dcomplex* a1;
	dcomplex* a2;
	dcomplex* g11;
	integer       j, g, k;
	integer       nG, nG_app;
	integer       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_opz_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = 0; j < k_minus_1; ++j )
	{
		nG_app = j + 1;

		for ( k = 0, g = nG_app - 1; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, rs_A,
			                     a2, rs_A );
		}
	}

	return FLA_SUCCESS;
}

