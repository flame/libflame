
#include "FLAME.h"

FLA_Error FLA_Apply_G_rf_asm_var3( FLA_Obj G, FLA_Obj A )
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

			FLA_Apply_G_rf_ass_var3( k_G,
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

			FLA_Apply_G_rf_asd_var3( k_G,
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

			FLA_Apply_G_rf_asc_var3( k_G,
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

			FLA_Apply_G_rf_asz_var3( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ass_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asd_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma23_k1;
	double    sigma23_k1;
	double    gamma34_k1;
	double    sigma34_k1;
	double    gamma12_k2;
	double    sigma12_k2;
	double    gamma23_k2;
	double    sigma23_k2;
	double*   a1;
	double*   a2;
	double*   a3;
	double*   a4;
	dcomplex* g23_k1;
	dcomplex* g34_k1;
	dcomplex* g12_k2;
	dcomplex* g23_k2;
	int       i, j, g, k;
	int       nG, nG_app;
	int       n_iter;
	int       n_left;
	int       k_minus_1;
	int       n_fuse;
	int       k_fuse;
	int       is_ident23_k1, is_ident34_k1;
	int       is_ident12_k2, is_ident23_k2;
	int       has_ident;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;
	n_fuse    = 2;
	k_fuse    = 2;

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
		nG_app = j + 2;
		n_iter = nG_app / k_fuse;
		//n_iter = nG_app % k_fuse;
		n_left = 1;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 2,
					                     a3, 2 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 2,
					                     a4, 2 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 2,
					                     a2, 2 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 2,
					                     a3, 2 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asd( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 2,
				                      a2, 2,
				                      a3, 2,
				                      a4, 2 );
			}
		}

		if ( n_left == 1 )
		{
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident34_k1 )
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 2,
				                     a4, 2 );
		}
	}

	// Pipeline stage

	for ( ; j < nG - 1; j += n_fuse )
	{
		nG_app = k_G;
		n_iter = nG_app / k_fuse;
		n_left = nG_app % k_fuse;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 2,
					                     a3, 2 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 2,
					                     a4, 2 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 2,
					                     a2, 2 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 2,
					                     a3, 2 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asd( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 2,
				                      a2, 2,
				                      a3, 2,
				                      a4, 2 );
			}
		}

		if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident23_k1 && is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 2,
				                     a3, 2 );
			}
			else if ( is_ident23_k1 && !is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 2,
				                     a4, 2 );
			}
			else
			{
				MAC_Apply_G_mx3_asd( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a2, 2,
				                     a3, 2,
				                     a4, 2 );
			}
		}
	}

	// Shutdown stage

	for ( j = nG % n_fuse; j < k_G; j += n_fuse )
	{
		g = nG - 1;
		k = j;

		//n_left = 1;
		//if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );

			if ( !is_ident23_k1 )
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 2,
				                     a3, 2 );
			++k;
			--g;
		}

		nG_app = k_minus_1 - j;
		n_iter = nG_app / k_fuse;
		n_left = nG_app % k_fuse;

		for ( i = 0; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 2,
					                     a3, 2 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 2,
					                     a4, 2 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 2,
					                     a2, 2 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asd( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 2,
					                     a3, 2 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asd( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 2,
				                      a2, 2,
				                      a3, 2,
				                      a4, 2 );
			}
		}

		if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident23_k1 && is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 2,
				                     a3, 2 );
			}
			else if ( is_ident23_k1 && !is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asd( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 2,
				                     a4, 2 );
			}
			else
			{
				MAC_Apply_G_mx3_asd( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a2, 2,
				                     a3, 2,
				                     a4, 2 );
			}
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A )
{
	double    one  = bl1_d1();
	double    zero = bl1_d0();
	double    gamma23_k1;
	double    sigma23_k1;
	double    gamma34_k1;
	double    sigma34_k1;
	double    gamma12_k2;
	double    sigma12_k2;
	double    gamma23_k2;
	double    sigma23_k2;
	dcomplex* a1;
	dcomplex* a2;
	dcomplex* a3;
	dcomplex* a4;
	dcomplex* g23_k1;
	dcomplex* g34_k1;
	dcomplex* g12_k2;
	dcomplex* g23_k2;
	int       i, j, g, k;
	int       nG, nG_app;
	int       n_iter;
	int       n_left;
	int       k_minus_1;
	int       n_fuse;
	int       k_fuse;
	int       is_ident23_k1, is_ident34_k1;
	int       is_ident12_k2, is_ident23_k2;
	int       has_ident;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;
	n_fuse    = 2;
	k_fuse    = 2;

	// Use the simple variant for nG < (k - 1) or k == 1.
	if ( nG < 2*k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_asz_var6( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = -1; j < k_minus_1; j += n_fuse )
	{
		nG_app = j + 2;
		n_iter = nG_app / k_fuse;
		//n_iter = nG_app % k_fuse;
		n_left = 1;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 1,
					                     a3, 1 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 1,
					                     a4, 1 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 1,
					                     a2, 1 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 1,
					                     a3, 1 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asz( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 1,
				                      a2, 1,
				                      a3, 1,
				                      a4, 1 );
			}
		}

		if ( n_left == 1 )
		{
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident34_k1 )
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 1,
				                     a4, 1 );
		}
	}

	// Pipeline stage

	for ( ; j < nG - 1; j += n_fuse )
	{
		nG_app = k_G;
		n_iter = nG_app / k_fuse;
		n_left = nG_app % k_fuse;

		for ( i = 0, k = 0, g = j; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 1,
					                     a3, 1 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 1,
					                     a4, 1 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 1,
					                     a2, 1 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 1,
					                     a3, 1 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asz( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 1,
				                      a2, 1,
				                      a3, 1,
				                      a4, 1 );
			}
		}

		if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident23_k1 && is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 1,
				                     a3, 1 );
			}
			else if ( is_ident23_k1 && !is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 1,
				                     a4, 1 );
			}
			else
			{
				MAC_Apply_G_mx3_asz( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a2, 1,
				                     a3, 1,
				                     a4, 1 );
			}
		}
	}

	// Shutdown stage

	for ( j = nG % n_fuse; j < k_G; j += n_fuse )
	{
		g = nG - 1;
		k = j;

		//n_left = 1;
		//if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );

			if ( !is_ident23_k1 )
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 1,
				                     a3, 1 );
			++k;
			--g;
		}

		nG_app = k_minus_1 - j;
		n_iter = nG_app / k_fuse;
		n_left = nG_app % k_fuse;

		for ( i = 0; i < n_iter; ++i, k += k_fuse, g -= n_fuse )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			g12_k2 = buff_G + (g - 1)*rs_G + (k + 1)*cs_G;
			g23_k2 = buff_G + (g    )*rs_G + (k + 1)*cs_G;
			a1     = buff_A + (g - 1)*cs_A;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;
			gamma12_k2 = g12_k2->real;
			sigma12_k2 = g12_k2->imag;
			gamma23_k2 = g23_k2->real;
			sigma23_k2 = g23_k2->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );
			is_ident12_k2 = ( gamma12_k2 == one && sigma12_k2 == zero );
			is_ident23_k2 = ( gamma23_k2 == one && sigma23_k2 == zero );
			has_ident     = ( is_ident23_k1 || is_ident34_k1 ||
			                  is_ident12_k2 || is_ident23_k2 );

			if      ( has_ident )
			{
				// Apply to pairs of columns as needed.

				if ( !is_ident23_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k1,
					                     &sigma23_k1,
					                     a2, 1,
					                     a3, 1 );

				if ( !is_ident34_k1 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma34_k1,
					                     &sigma34_k1,
					                     a3, 1,
					                     a4, 1 );

				if ( !is_ident12_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma12_k2,
					                     &sigma12_k2,
					                     a1, 1,
					                     a2, 1 );

				if ( !is_ident23_k2 )
					MAC_Apply_G_mx2_asz( m_A,
					                     &gamma23_k2,
					                     &sigma23_k2,
					                     a2, 1,
					                     a3, 1 );
			}
			else
			{
				// Apply to all four columns.

				MAC_Apply_G_mx4s_asz( m_A,
				                      &gamma23_k1,
				                      &sigma23_k1,
				                      &gamma34_k1,
				                      &sigma34_k1,
				                      &gamma12_k2,
				                      &sigma12_k2,
				                      &gamma23_k2,
				                      &sigma23_k2,
				                      a1, 1,
				                      a2, 1,
				                      a3, 1,
				                      a4, 1 );
			}
		}

		if ( n_left == 1 )
		{
			g23_k1 = buff_G + (g    )*rs_G + (k    )*cs_G;
			g34_k1 = buff_G + (g + 1)*rs_G + (k    )*cs_G;
			a2     = buff_A + (g    )*cs_A;
			a3     = buff_A + (g + 1)*cs_A;
			a4     = buff_A + (g + 2)*cs_A;

			gamma23_k1 = g23_k1->real;
			sigma23_k1 = g23_k1->imag;
			gamma34_k1 = g34_k1->real;
			sigma34_k1 = g34_k1->imag;

			is_ident23_k1 = ( gamma23_k1 == one && sigma23_k1 == zero );
			is_ident34_k1 = ( gamma34_k1 == one && sigma34_k1 == zero );

			if ( !is_ident23_k1 && is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     a2, 1,
				                     a3, 1 );
			}
			else if ( is_ident23_k1 && !is_ident34_k1 )
			{
				MAC_Apply_G_mx2_asz( m_A,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a3, 1,
				                     a4, 1 );
			}
			else
			{
				MAC_Apply_G_mx3_asz( m_A,
				                     &gamma23_k1,
				                     &sigma23_k1,
				                     &gamma34_k1,
				                     &sigma34_k1,
				                     a2, 1,
				                     a3, 1,
				                     a4, 1 );
			}
		}
	}

	return FLA_SUCCESS;
}

