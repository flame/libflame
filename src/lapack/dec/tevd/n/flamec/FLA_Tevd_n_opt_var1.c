/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_n_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Datatype datatype;
	int          m_A, m_U, n_G;
	int          inc_d;
	int          inc_e;
	int          rs_G, cs_G;

	datatype = FLA_Obj_datatype( U );

	m_A       = FLA_Obj_vector_dim( d );
	m_U       = FLA_Obj_vector_dim( d );
	n_G       = FLA_Obj_width( G );

	inc_d     = FLA_Obj_vector_inc( d );
	inc_e     = FLA_Obj_vector_inc( e );
	
	rs_G      = FLA_Obj_row_stride( G );
	cs_G      = FLA_Obj_col_stride( G );

/*
FLA_Obj de, deL, deR, deLT, deLB;
FLA_Obj_create( FLA_DOUBLE, m_A, 2, 0, 0, &de );
FLA_Part_1x2( de,  &deL, &deR,   1, FLA_LEFT );
FLA_Part_2x1( deL, &deLT,
                   &deLB,   1, FLA_BOTTOM );
FLA_Copy( d, deR );
FLA_Copy( e, deLT );
FLA_Set( FLA_ZERO, deLB );
//FLA_Obj_show( "de", de, "%21.17e", "" );
*/

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );

			r_val = FLA_Tevd_n_ops_var1( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );

			r_val = FLA_Tevd_n_opd_var1( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );

			r_val = FLA_Tevd_n_opc_var1( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );

			r_val = FLA_Tevd_n_opz_var1( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G );

			break;
		}
	}
/*
FLA_Copy( d, deR );
FLA_Copy( e, deLT );
FLA_Set( FLA_ZERO, deLB );
FLA_Sort( FLA_FORWARD, deR );
FLA_Obj_show( "de after", de, "%21.17e", "" );
double eps   = FLA_Mach_params_opd( FLA_MACH_EPS );
printf( "epsilon = %21.17e\n", eps );
FLA_Obj_free( &de );
*/
	return r_val;
}



FLA_Error FLA_Tevd_n_ops_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G )
{
	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_n_opd_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_n_opc_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_n_opz_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G )
{
	dcomplex  one  = bl1_z1();
	double    rone = bl1_d1();

	double    eps;
	double    eps2;
	double    safmin;
	double    ssfmin;
	double    safmax;
	double    ssfmax;

	dcomplex* G;
	double*   d1;
	double*   e1;
	int       r_val;
	int       done;
	int       m_G_sweep_max;
	int       ij_begin;
	int       ijTL, ijBR;
	int       m_A11;
	int       n_iter_perf;
	int       total_deflations;
	int       n_deflations;
	int       n_iter_prev;
	int       n_iter_perf_sweep_max;

	// Initialize some numerical constants.
	eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	eps2   = FLA_Mach_params_opd( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );
	safmax = rone / safmin;
	ssfmax = sqrt( safmax ) / 3.0;
	ssfmin = sqrt( safmin ) / eps2;

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Iterate until the matrix has completely deflated.
	for ( total_deflations = 0; done != TRUE; )
	{

		// Initialize G to contain only identity rotations.
		bl1_zsetm( m_G_sweep_max,
		           n_G,
		           &one,
		           buff_G, rs_G, cs_G );

		// Keep track of the maximum number of iterations performed in the
		// current sweep. This is used when applying the sweep's Givens
		// rotations.
		n_iter_perf_sweep_max = 0;

		// Perform a sweep: Move through the matrix and perform a tridiagonal
		// EVD on each non-zero submatrix that is encountered. During the
		// first time through, ijTL will be 0 and ijBR will be m_A - 1.
		for ( ij_begin = 0; ij_begin < m_A; )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_n_opz_var1: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first submatrix along the diagonal that is
			// bounded by zeroes (or endpoints of the matrix). If no
			// submatrix is found (ie: if the entire subdiagonal is zero
			// then FLA_FAILURE is returned. This function also inspects
			// subdiagonal elements for proximity to zero. If a given
			// element is close enough to zero, then it is deemed
			// converged and manually set to zero.
			r_val = FLA_Tevd_find_submatrix_opd( m_A,
			                                     ij_begin,
			                                     buff_d, inc_d,
			                                     buff_e, inc_e,
			                                     &ijTL,
			                                     &ijBR );

			// Verify that a submatrix was found. If one was not found,
			// then we are done with the current sweep. Furthermore, if
			// a submatrix was not found AND we began our search at the
			// beginning of the matrix (ie: ij_begin == 0), then the
			// matrix has completely deflated and so we are done with
			// Francis step iteration.
			if ( r_val == FLA_FAILURE )
			{
				if ( ij_begin == 0 )
				{
#ifdef PRINTF
printf( "FLA_Tevd_n_opz_var1: subdiagonal is completely zero.\n" );
printf( "FLA_Tevd_n_opz_var1: Francis iteration is done!\n" );
#endif
					done = TRUE;
				}

				// Break out of the current sweep so we can apply the last
				// remaining Givens rotations.
				break;
			}

			// If we got this far, then:
			//   (a) ijTL refers to the index of the first non-zero
			//       subdiagonal along the diagonal, and
			//   (b) ijBR refers to either:
			//       - the first zero element that occurs after ijTL, or
			//       - the the last diagonal element.
			// Note that ijTL and ijBR also correspond to the first and
			// last diagonal elements of the submatrix of interest. Thus,
			// we may compute the dimension of this submatrix as:
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Tevd_n_opz_var1: ij_begin = %d\n", ij_begin );
printf( "FLA_Tevd_n_opz_var1: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_n_opz_var1: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_n_opz_var1: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next submatrix
			// search in the current sweep.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;

			// Compute the 1-norm (which equals the infinity norm since the
			// matrix is tridiagonal and symmetric) and perform appropriate
			// scaling if necessary.
/*
			FLA_Norm1_tridiag( m_A11,
			                   d1, inc_d,
			                   e1, inc_e,
			                   &norm );
*/	

			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever a split occurs. Iteration continues
			// as long as:
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Tevd_iteracc_n_opd_var1( m_A11,
			                                            n_G,
			                                            ijTL,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            &n_iter_perf );

			// Record the number of deflations that were observed.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Tevd_n_opz_var1: deflations observed       = %d\n", n_deflations );
printf( "FLA_Tevd_n_opz_var1: total deflations observed = %d\n", total_deflations );
printf( "FLA_Tevd_n_opz_var1: num iterations performed  = %d\n", n_iter_perf );
#endif

			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G we can apply and safely include all
			// non-identity rotations that were computed during the
			// eigenvalue searches.
			m_G_sweep_max = ijBR;

			// Make sure we haven't exceeded our maximum iteration count.
			if ( n_iter_prev >= m_A * n_iter_max )
			{
#ifdef PRINTF
printf( "FLA_Tevd_n_opz_var1: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
				FLA_Abort();
				//return FLA_FAILURE;
			}
		}

		// The sweep is complete.

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
printf( "FLA_Tevd_n_opz_var1: total number of iterations performed: %d\n", n_iter_prev );
#endif
	}

	//return FLA_SUCCESS;
	return n_iter_prev;
}

