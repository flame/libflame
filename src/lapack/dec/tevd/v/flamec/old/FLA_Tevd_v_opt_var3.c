/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_v_opt_var3( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U, dim_t b_alg )
{
	FLA_Datatype datatype;
	integer          m_A, m_U, n_G;
	integer          n_G_extra;
	integer          inc_d;
	integer          inc_e;
	integer          rs_G, cs_G;
	integer          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A       = FLA_Obj_vector_dim( d );
	m_U       = FLA_Obj_length( U );
	n_G       = FLA_Obj_width( G ) - n_iter_max;
	n_G_extra = n_iter_max;

	inc_d     = FLA_Obj_vector_inc( d );
	inc_e     = FLA_Obj_vector_inc( e );
	
	rs_G      = FLA_Obj_row_stride( G );
	cs_G      = FLA_Obj_col_stride( G );

	rs_U      = FLA_Obj_row_stride( U );
	cs_U      = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_U = FLA_FLOAT_PTR( U );

			FLA_Tevd_v_ops_var3( m_A,
			                     m_U,
			                     n_G,
			                     n_G_extra,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U,
			                     b_alg );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			FLA_Tevd_v_opd_var3( m_A,
			                     m_U,
			                     n_G,
			                     n_G_extra,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U,
			                     b_alg );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			FLA_Tevd_v_opc_var3( m_A,
			                     m_U,
			                     n_G,
			                     n_G_extra,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U,
			                     b_alg );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			FLA_Tevd_v_opz_var3( m_A,
			                     m_U,
			                     n_G,
			                     n_G_extra,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U,
			                     b_alg );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var3( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_opd_var3( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opc_var3( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opz_var3( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	dcomplex  one        = bl1_z1();
	double    zero       = bl1_d0();
	integer       n_iter_max = n_G_extra;

	dcomplex* G;
	double*   d1;
	double*   e1;
	integer       total_eigs_found;
	integer       n_eigs_found;
	integer       n_iter_perf;
	integer       n_U_apply;
	integer       m_G_sweep_max;
	integer       ij_begin;
	integer       ijTL, ijBR;
	integer       m_A11;

	// Initialze a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Iterate until we've found all eigenvalues.
	for ( total_eigs_found = 0; total_eigs_found < m_A; )
	{

		// Initialize G to contain only identity rotations.
		bl1_zsetm( m_G_sweep_max,
		           n_G + n_iter_max,
		           &one,
		           buff_G, rs_G, cs_G );

		// Perform a recursive tridiagonal decomposition on each submatrix
		// that falls between zeroes (or endpoints) along the subdiagonal.
		// During the first time through, ijTL will be 0 and ijBR will be
		// m_A - 1.
		for ( ij_begin = 0; ij_begin < m_A;  )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_v_opz_var3: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first non-zero subdiagonal element.
			for ( ijTL = ij_begin; ijTL < m_A; ++ijTL )
			{
				if ( buff_e[ (ijTL)*inc_e ] != zero )
				{
#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: found nonzero e[%d] = %21.19f\n", ijTL, buff_e[ (ijTL)*inc_e ] );
#endif
					break;
				}
			}
			// If we didn't find any non-zeros, then we are done with this
			// sweep of subproblems.
			if ( ijTL == m_A )
			{
#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: done with current sweep.\n" );
#endif
				break;
			}

			// Search for the first zero subdiagonal element after ijTL.
			// Once we find this element (ijBR) of the subdiagonal, this
			// corresponds to the adjacent element of d, d[ijBR], that is
			// non-zero.
			for ( ijBR = ijTL; ijBR < m_A - 1; ++ijBR )
			{
				if ( buff_e[ (ijBR)*inc_e ] == zero )
				{
#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: found zero e[%d] = %21.19f\n", ijBR, buff_e[ (ijBR)*inc_e ] );
#endif
					break;
				}
			}

			// At this point, ijTL and ijBR refer to the indices of the first
			// and last elements, respectively, of the current tridiagonal
			// submatrix.

			// Determine the dimension of the submatrix that was just
			// partitioned.
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Tevd_v_opt_var3: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_v_opt_var3: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_v_opt_var3: m_A11    = %d\n", m_A11 );
printf( "FLA_Tevd_v_opt_var3: ij_begin = %d\n", ij_begin );
#endif

			// Adjust ij_begin, which gets us ready for the next subproblem, if
			// there is one.
			ij_begin = ijBR + 1;

		
			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever possible. A new eigenvalue search is
			// performed as long as
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;

			n_eigs_found = FLA_Tevd_v_opd_var3r( m_A11,
			                                     n_G,
			                                     n_iter_max,
			                                     ijTL,
			                                     d1, inc_d,
			                                     e1, inc_e,
			                                     G,  rs_G, cs_G,
			                                     &n_iter_perf );

			// Record the number of eigenvalues found.
			total_eigs_found += n_eigs_found;

#ifdef PRINTF
printf( "FLA_Tevd_v_opt_var3: num eigs found       = %d\n", n_eigs_found );
printf( "FLA_Tevd_v_opt_var3: total num eigs found = %d\n", total_eigs_found );
printf( "FLA_Tevd_v_opt_var3: num iterations       = %d\n", n_iter_perf );
#endif


#ifdef PRINTF
printf( "FLA_Tevd_v_opt_var3: applying %d sets of Givens rotations\n", n_iter_perf );
printf( "FLA_Tevd_v_opt_var3: m_U = %d\n", m_U );
printf( "FLA_Tevd_v_opt_var3: m_A = %d\n", m_A );
#endif
			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G we can apply and safely include all
			// non-identity rotations that were computed during the
			// eigenvalue searches.
			m_G_sweep_max = ijBR;

			// Recall that the number of columns of U to which we apply
			// rotations is one more than the number of rotations.
			n_U_apply = m_G_sweep_max + 1;

			// Apply the Givens rotations that were computed as part of
			// the previous batch of iterations.
			FLA_Apply_G_rf_blz_var3( n_iter_perf,
			                         m_U,
			                         n_U_apply,
			                         buff_G, rs_G, cs_G,
			                         buff_U, rs_U, cs_U,
			                         b_alg );
		}
	}

	return FLA_SUCCESS;
}

