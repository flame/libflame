
#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_opt_var3( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of G:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  This variant iterates in pipelined, non-overlapping fashion and
  applies rotations to two columns at a time.

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

			FLA_Apply_G_rf_ops_var3( k_G,
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

			FLA_Apply_G_rf_opd_var3( k_G,
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

			FLA_Apply_G_rf_opc_var3( k_G,
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

			FLA_Apply_G_rf_opz_var3( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ops_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opd_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma;
	double    sigma;
	double*   a1;
	double*   a2;
	dcomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       twok_minus_1;
	int       offset;

	twok_minus_1 = 2 * ( k_G - 1 );
	nG           = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < twok_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_opd_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}

	// Start-up phase.

	for ( j = 0; j < twok_minus_1; ++j )
	{
		nG_app = j / 2 + 1;
		offset = j % 2;

		for ( g = 0, k = nG_app - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (offset + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (offset + 2*g    )*cs_A;
			a2    = buff_A + (offset + 2*g + 1)*cs_A;

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

	for ( j = 0; j < nG - twok_minus_1; ++j )
	{
		nG_app = k_G;

		for ( g = 0, k = nG_app - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (j + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (j + 2*g    )*cs_A;
			a2    = buff_A + (j + 2*g + 1)*cs_A;

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

	for ( j = nG - twok_minus_1; j < nG; ++j )
	{
		nG_app = ( nG - j + 1 ) / 2;

		for ( g = 0, k = k_G - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (j + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (j + 2*g    )*cs_A;
			a2    = buff_A + (j + 2*g + 1)*cs_A;

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

FLA_Error FLA_Apply_G_rf_opc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_opz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma;
	double    sigma;
	dcomplex* a1;
	dcomplex* a2;
	dcomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       twok_minus_1;
	int       offset;

	twok_minus_1 = 2 * ( k_G - 1 );
	nG           = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < twok_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_opz_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}

	// Start-up phase.

	for ( j = 0; j < twok_minus_1; ++j )
	{
		nG_app = j / 2 + 1;
		offset = j % 2;

		for ( g = 0, k = nG_app - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (offset + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (offset + 2*g    )*cs_A;
			a2    = buff_A + (offset + 2*g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Pipeline stage

	for ( j = 0; j < nG - twok_minus_1; ++j )
	{
		nG_app = k_G;

		for ( g = 0, k = nG_app - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (j + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (j + 2*g    )*cs_A;
			a2    = buff_A + (j + 2*g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Shutdown stage

	for ( j = nG - twok_minus_1; j < nG; ++j )
	{
		nG_app = ( nG - j + 1 ) / 2;

		for ( g = 0, k = k_G - 1; g < nG_app; ++g, --k )
		{
			g11   = buff_G + (j + 2*g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (j + 2*g    )*cs_A;
			a2    = buff_A + (j + 2*g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_opz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	return FLA_SUCCESS;
}

