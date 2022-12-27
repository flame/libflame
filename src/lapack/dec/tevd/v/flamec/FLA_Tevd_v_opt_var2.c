/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, dim_t b_alg )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Datatype datatype;
	integer          m_A, m_U, n_G;
	integer          inc_d;
	integer          inc_e;
	integer          rs_G, cs_G;
	integer          rs_R, cs_R;
	integer          rs_U, cs_U;
	integer          rs_W, cs_W;

	datatype = FLA_Obj_datatype( U );

	m_A       = FLA_Obj_vector_dim( d );
	m_U       = FLA_Obj_length( U );
	n_G       = FLA_Obj_width( G );

	inc_d     = FLA_Obj_vector_inc( d );
	inc_e     = FLA_Obj_vector_inc( e );
	
	rs_G      = FLA_Obj_row_stride( G );
	cs_G      = FLA_Obj_col_stride( G );

	rs_R      = FLA_Obj_row_stride( R );
	cs_R      = FLA_Obj_col_stride( R );

	rs_W      = FLA_Obj_row_stride( W );
	cs_W      = FLA_Obj_col_stride( W );

	rs_U      = FLA_Obj_row_stride( U );
	cs_U      = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_R = FLA_FLOAT_PTR( R );
			float*    buff_W = FLA_FLOAT_PTR( W );
			float*    buff_U = FLA_FLOAT_PTR( U );

			r_val = FLA_Tevd_v_ops_var2( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_R, rs_R, cs_R,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_R = FLA_DOUBLE_PTR( R );
			double*   buff_W = FLA_DOUBLE_PTR( W );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			r_val = FLA_Tevd_v_opd_var2( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_R, rs_R, cs_R,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_R = FLA_FLOAT_PTR( R );
			scomplex* buff_W = FLA_COMPLEX_PTR( W );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			r_val = FLA_Tevd_v_opc_var2( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_R, rs_R, cs_R,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_R = FLA_DOUBLE_PTR( R );
			dcomplex* buff_W = FLA_DOUBLE_COMPLEX_PTR( W );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			r_val = FLA_Tevd_v_opz_var2( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_R, rs_R, cs_R,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}
	}

	return r_val;
}



FLA_Error FLA_Tevd_v_ops_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_R, integer rs_R, integer cs_R,
                               float*    buff_W, integer rs_W, integer cs_W,
                               float*    buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_opd_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_R, integer rs_R, integer cs_R,
                               double*   buff_W, integer rs_W, integer cs_W,
                               double*   buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	dcomplex  one   = bl1_z1();
	double    rone  = bl1_d1();
	double    rzero = bl1_d0();

	dcomplex* G;
	double*   d1;
	double*   e1;
	integer       r_val;
	integer       done;
	integer       m_G_sweep_max;
	integer       ij_begin;
	integer       ijTL, ijBR;
	integer       m_A11;
	integer       n_iter_perf;
	integer       n_U_apply;
	integer       total_deflations;
	integer       n_deflations;
	integer       n_iter_prev;
	integer       n_iter_perf_sweep_max;

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Initialize R to identity.
	bl1_dident( m_A,
	            buff_R, rs_R, cs_R );

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
		for ( ij_begin = 0; ij_begin < m_A;  )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_v_opd_var2: beginning new sweep (ij_begin = %d)\n", ij_begin );
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
printf( "FLA_Tevd_v_opd_var2: subdiagonal is completely zero.\n" );
printf( "FLA_Tevd_v_opd_var2: Francis iteration is done!\n" );
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
printf( "FLA_Tevd_v_opd_var2: ij_begin = %d\n", ij_begin );
printf( "FLA_Tevd_v_opd_var2: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_v_opd_var2: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_v_opd_var2: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next subproblem, if
			// there is one.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;

			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever a split occurs. Iteration continues
			// as long as
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Tevd_iteracc_v_opd_var1( m_A11,
			                                            n_G,
			                                            ijTL,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            G,  rs_G, cs_G,
			                                            &n_iter_perf );

			// Record the number of deflations that we observed.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = fla_max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var2: deflations observed       = %d\n", n_deflations );
printf( "FLA_Tevd_v_opd_var2: total deflations observed = %d\n", total_deflations );
printf( "FLA_Tevd_v_opd_var2: num iterations            = %d\n", n_iter_perf );
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
printf( "FLA_Tevd_v_opd_var2: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
				FLA_Abort();
				//return FLA_FAILURE;
			}
		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.


		// Recall that the number of columns of U to which we apply
		// rotations is one more than the number of rotations.
		n_U_apply = m_G_sweep_max + 1;

		// Apply the Givens rotations that were computed as part of
		// the previous batch of iterations.
		//FLA_Apply_G_rf_bld_var8b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          m_U,
		                          n_U_apply,
		                          n_iter_prev,
		                          buff_G, rs_G, cs_G,
		                          buff_R, rs_R, cs_R,
		                          b_alg );

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var2: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
#endif

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;
	}

	// Copy the contents of Q to temporary storage.
	bl1_dcopymt( BLIS1_NO_TRANSPOSE,
	             m_A,
	             m_A,
	             buff_U, rs_U, cs_U,
	             buff_W, rs_W, cs_W );


	// Multiply Q by R, overwriting U.
	bl1_dgemm( BLIS1_NO_TRANSPOSE,
	           BLIS1_NO_TRANSPOSE,
	           m_A,
	           m_A,
	           m_A,
	           &rone,
	           ( double* )buff_W, rs_W, cs_W,
	                      buff_R, rs_R, cs_R,
	           &rzero,
	           ( double* )buff_U, rs_U, cs_U );

	return n_iter_prev;
}

FLA_Error FLA_Tevd_v_opc_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_R, integer rs_R, integer cs_R,
                               scomplex* buff_W, integer rs_W, integer cs_W,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opz_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_R, integer rs_R, integer cs_R,
                               dcomplex* buff_W, integer rs_W, integer cs_W,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg )
{
	dcomplex  one   = bl1_z1();
	double    rone  = bl1_d1();
	double    rzero = bl1_d0();

	dcomplex* G;
	double*   d1;
	double*   e1;
	integer       r_val;
	integer       done;
	integer       m_G_sweep_max;
	integer       ij_begin;
	integer       ijTL, ijBR;
	integer       m_A11;
	integer       n_iter_perf;
	integer       n_U_apply;
	integer       total_deflations;
	integer       n_deflations;
	integer       n_iter_prev;
	integer       n_iter_perf_sweep_max;

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Initialize R to identity.
	bl1_dident( m_A,
	            buff_R, rs_R, cs_R );

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
		for ( ij_begin = 0; ij_begin < m_A;  )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_v_opz_var2: beginning new sweep (ij_begin = %d)\n", ij_begin );
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
printf( "FLA_Tevd_v_opz_var2: subdiagonal is completely zero.\n" );
printf( "FLA_Tevd_v_opz_var2: Francis iteration is done!\n" );
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
printf( "FLA_Tevd_v_opz_var2: ij_begin = %d\n", ij_begin );
printf( "FLA_Tevd_v_opz_var2: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_v_opz_var2: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_v_opz_var2: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next subproblem, if
			// there is one.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;

			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever a split occurs. Iteration continues
			// as long as
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Tevd_iteracc_v_opd_var1( m_A11,
			                                            n_G,
			                                            ijTL,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            G,  rs_G, cs_G,
			                                            &n_iter_perf );

			// Record the number of deflations that we observed.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = fla_max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var2: deflations observed       = %d\n", n_deflations );
printf( "FLA_Tevd_v_opz_var2: total deflations observed = %d\n", total_deflations );
printf( "FLA_Tevd_v_opz_var2: num iterations            = %d\n", n_iter_perf );
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
printf( "FLA_Tevd_v_opz_var2: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
				FLA_Abort();
				//return FLA_FAILURE;
			}
		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.


		// Recall that the number of columns of U to which we apply
		// rotations is one more than the number of rotations.
		n_U_apply = m_G_sweep_max + 1;

		// Apply the Givens rotations that were computed as part of
		// the previous batch of iterations.
		//FLA_Apply_G_rf_bld_var8b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          m_U,
		                          n_U_apply,
		                          n_iter_prev,
		                          buff_G, rs_G, cs_G,
		                          buff_R, rs_R, cs_R,
		                          b_alg );

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var2: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
#endif

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;
	}

	// Copy the contents of Q to temporary storage.
	bl1_zcopymt( BLIS1_NO_TRANSPOSE,
	             m_A,
	             m_A,
	             buff_U, rs_U, cs_U,
	             buff_W, rs_W, cs_W );


	// Multiply Q by R, overwriting U.
	bl1_dgemm( BLIS1_NO_TRANSPOSE,
	           BLIS1_NO_TRANSPOSE,
	           2*m_A,
	           m_A,
	           m_A,
	           &rone,
	           ( double* )buff_W, rs_W, 2*cs_W,
	                      buff_R, rs_R,   cs_R,
	           &rzero,
	           ( double* )buff_U, rs_U, 2*cs_U );

	return n_iter_prev;
}

