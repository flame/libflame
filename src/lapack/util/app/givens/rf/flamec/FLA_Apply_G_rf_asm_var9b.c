/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_asm_var9b( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of C and S:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

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

			FLA_Apply_G_rf_ass_var9b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

			FLA_Apply_G_rf_asd_var9b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asc_var9b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asz_var9b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ass_var9b( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   integer       i_k,
                                   integer       iTL,
                                   scomplex* buff_G, integer rs_G, integer cs_G,
                                   float*    buff_A, integer rs_A, integer cs_A )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asd_var9b( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   integer       i_k,
                                   integer       iTL,
                                   dcomplex* buff_G, integer rs_G, integer cs_G,
                                   double*   buff_A, integer rs_A, integer cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma12;
	double    sigma12;
	double    gamma23;
	double    sigma23;
	double*   a1;
	double*   a2;
	double*   a3;
	dcomplex* g12;
	dcomplex* g23;
	integer       i, j, g, k;
	integer       nG, nG_app;
	integer       n_iter;
	integer       n_left;
	integer       k_minus_1;
	integer       n_fuse;
	integer       is_ident12, is_ident23;
	integer       m_app;


	k_minus_1 = k_G - 1;
	nG        = n_A - 1;
	n_fuse    = 2;

	// Use the simple variant for nG < (k - 1) or k == 1.
	if ( nG < 2*k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_asd_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = -1; j < k_minus_1; j += n_fuse )
	{
		nG_app = j + 1;
		n_iter = nG_app;
		n_left = 1;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, ++k, --g )
		{
			g12   = buff_G + (g    )*rs_G + (k    )*cs_G;
			g23   = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;
			a3    = buff_A + (g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

			m_app = min( i_k + 3 + j - iTL, m_A );
			m_app = max( m_app, 0 );

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 1,
				                     a2, 1 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 1,
				                     a3, 1 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 1,
				                     a2, 1,
				                     a3, 1 );
			}
		}

		if ( n_left == 1 )
		{
			g23   = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a2    = buff_A + (g + 1)*cs_A;
			a3    = buff_A + (g + 2)*cs_A;

			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident23 = ( gamma23 == one && sigma23 == zero );

			m_app = min( i_k + 3 + j - iTL, m_A );
			m_app = max( m_app, 0 );

			if ( !is_ident23 )
				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 1,
				                     a3, 1 );
		}
	}

	// Pipeline stage

	for ( ; j < nG - 1; j += n_fuse )
	{
		nG_app = k_G;
		n_iter = nG_app;
		n_left = 0;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, ++k, --g )
		{
			g12   = buff_G + (g    )*rs_G + (k    )*cs_G;
			g23   = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;
			a3    = buff_A + (g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

			m_app = min( i_k + 3 + j - iTL, m_A );
			m_app = max( m_app, 0 );

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 1,
				                     a2, 1 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 1,
				                     a3, 1 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 1,
				                     a2, 1,
				                     a3, 1 );
			}
		}
	}

	// Shutdown stage

	for ( j = nG % n_fuse; j < k_G; j += n_fuse )
	{
		g = nG - 1;
		k = j;

		n_left = 1;
		if ( n_left == 1 )
		{
			g12   = buff_G + (g    )*rs_G + (k    )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );

			m_app = m_A;

			if ( !is_ident12 )
				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 1,
				                     a2, 1 );
			++k;
			--g;
		}

		nG_app = k_minus_1 - j;
		n_iter = nG_app;

		for ( i = 0; i < n_iter; ++i, ++k, --g )
		{
			g12   = buff_G + (g    )*rs_G + (k    )*cs_G;
			g23   = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;
			a3    = buff_A + (g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

			m_app = m_A;

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 1,
				                     a2, 1 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 1,
				                     a3, 1 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 1,
				                     a2, 1,
				                     a3, 1 );
			}
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asc_var9b( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   integer       i_k,
                                   integer       iTL,
                                   scomplex* buff_G, integer rs_G, integer cs_G,
                                   scomplex* buff_A, integer rs_A, integer cs_A )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asz_var9b( integer       k_G,
                                   integer       m_A,
                                   integer       n_A,
                                   integer       i_k,
                                   integer       iTL,
                                   dcomplex* buff_G, integer rs_G, integer cs_G,
                                   dcomplex* buff_A, integer rs_A, integer cs_A )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

