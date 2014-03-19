/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_opt_var1( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of G:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  This variant iterates naively and applies rotations to two columns
  at a time.

  -FGVZ
*/
{
	FLA_Datatype datatype;
	int          k_G, m_A, n_A;
	int          rs_G, cs_G;
	int          rs_A, cs_A;

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

			FLA_Apply_G_rf_ops_var1( k_G,
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

			FLA_Apply_G_rf_opd_var1( k_G,
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

			FLA_Apply_G_rf_opc_var1( k_G,
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

			FLA_Apply_G_rf_opz_var1( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ops_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A )
{
	float     one    = bl1_s1();
	float     zero   = bl1_s0();
	int       nG_app = n_A - 1;
	int       l, j;
	float     gamma;
	float     sigma;
	float*    a1;
	float*    a2;
	scomplex* g1;
	scomplex* g11;

	g1 = buff_G;

	for ( l = 0; l < k_G; ++l )
	{
		a1 = buff_A;
		a2 = buff_A + cs_A;
		g11 = g1;

		for ( j = 0; j < nG_app; ++j )
		{
			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma != one || sigma != zero )
			{
				MAC_Apply_G_mx2_ops( m_A,
				                     &gamma,
				                     &sigma,
				                     a1, rs_A,
				                     a2, rs_A );
			}

			a1 += cs_A;
			a2 += cs_A;
			g11 += rs_G;
		}

		g1 += cs_G;
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opd_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A )
{
	double    one    = bl1_d1();
	double    zero   = bl1_d0();
	int       nG_app = n_A - 1;
	int       l, j;
	double    gamma;
	double    sigma;
	double*   a1;
	double*   a2;
	dcomplex* g1;
	dcomplex* g11;

	g1 = buff_G;

	for ( l = 0; l < k_G; ++l )
	{
		a1 = buff_A;
		a2 = buff_A + cs_A;
		g11 = g1;

		for ( j = 0; j < nG_app; ++j )
		{
			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma != one || sigma != zero )
			{
				MAC_Apply_G_mx2_opd( m_A,
				                     &gamma,
				                     &sigma,
				                     a1, rs_A,
				                     a2, rs_A );
			}

			a1 += cs_A;
			a2 += cs_A;
			g11 += rs_G;
		}

		g1 += cs_G;
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opc_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A )
{
	float     one    = bl1_s1();
	float     zero   = bl1_s0();
	int       nG_app = n_A - 1;
	int       l, j;
	float     gamma;
	float     sigma;
	scomplex* a1;
	scomplex* a2;
	scomplex* g1;
	scomplex* g11;

	g1 = buff_G;

	for ( l = 0; l < k_G; ++l )
	{
		a1 = buff_A;
		a2 = buff_A + cs_A;
		g11 = g1;

		for ( j = 0; j < nG_app; ++j )
		{
			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma != one || sigma != zero )
			{
				MAC_Apply_G_mx2_opc( m_A,
				                     &gamma,
				                     &sigma,
				                     a1, rs_A,
				                     a2, rs_A );
			}

			a1 += cs_A;
			a2 += cs_A;
			g11 += rs_G;
		}

		g1 += cs_G;
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_opz_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A )
{
	double    one    = bl1_d1();
	double    zero   = bl1_d0();
	int       nG_app = n_A - 1;
	int       l, j;
	double    gamma;
	double    sigma;
	dcomplex* a1;
	dcomplex* a2;
	dcomplex* g1;
	dcomplex* g11;

	g1 = buff_G;

	for ( l = 0; l < k_G; ++l )
	{
		a1 = buff_A;
		a2 = buff_A + cs_A;
		g11 = g1;

		for ( j = 0; j < nG_app; ++j )
		{
			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma != one || sigma != zero )
			{
				MAC_Apply_G_mx2_opz( m_A,
				                     &gamma,
				                     &sigma,
				                     a1, rs_A,
				                     a2, rs_A );
			}

			a1 += cs_A;
			a2 += cs_A;
			g11 += rs_G;
		}

		g1 += cs_G;
	}

	return FLA_SUCCESS;
}

