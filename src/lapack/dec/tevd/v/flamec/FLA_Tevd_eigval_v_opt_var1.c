/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_eigval_v_opt_var1( FLA_Obj G, FLA_Obj d, FLA_Obj e, FLA_Obj k )
{
	FLA_Datatype datatype;
	int          m_A, n_G;
	int          rs_G, cs_G;
	int          inc_d;
	int          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );
	n_G      = FLA_Obj_width( G );

	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			int*      buff_k = FLA_INT_PTR( k );

			FLA_Tevd_eigval_v_ops_var1( m_A,
			                            n_G,
			                            buff_G, rs_G, cs_G,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_k );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			int*      buff_k = FLA_INT_PTR( k );

			FLA_Tevd_eigval_v_opd_var1( m_A,
			                            n_G,
			                            buff_G, rs_G, cs_G,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_k );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_eigval_v_ops_var1( int       m_A,
                                      int       n_G,
                                      scomplex* buff_G, int rs_G, int cs_G,
                                      float*    buff_d, int inc_d, 
                                      float*    buff_e, int inc_e,
                                      int*      n_iter )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_eigval_v_opd_var1( int       m_A,
                                      int       n_G,
                                      dcomplex* buff_G, int rs_G, int cs_G,
                                      double*   buff_d, int inc_d, 
                                      double*   buff_e, int inc_e,
                                      int*      n_iter )
{
	FLA_Error r_val;
	double    eps;
	double    safmin;
	double*   e_last;
	double*   d_last;
	double*   d_last_m1;
	double    shift;
	int       k;
	int       n_iter_allowed = n_G;

	// Query epsilon and safmin, which are used in the test for convergence.
	eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );

	// Initialize a pointer to the last sub-diagonal element and two
	// more to the last and second last
	e_last    = &buff_e[ (m_A-2)*inc_e ];
	d_last_m1 = &buff_d[ (m_A-2)*inc_d ];
	d_last    = &buff_d[ (m_A-1)*inc_d ];

	for ( k = 0; k < n_iter_allowed; ++k )
	{
		dcomplex* g1 = buff_G + (k  )*cs_G;

		/*------------------------------------------------------------*/

		// If we've converged, record k and return index of eigenvalue found.
		// The reason we check before the Francis step (rather than after)
		// is so we correctly handle situations where the last diagonal
		// element has already converged from previous eigenvalue searches
		// and thus no iteration is necessary. If we checked after the
		// Francis step, we would have unnecessarily executed an additional
		// Francis step's worth of rotations with a sub-optimal shift (since
		// it would be using a 2x2 that was not "centered" properly).
		if ( MAC_Tevd_eigval_converged_opd( eps, safmin, *d_last_m1, *e_last, *d_last ) )
		{
			*e_last = 0.0;
			*n_iter = k;
			return m_A - 1;
		}

		// Compute a Wilkinson shift with the last 2x2 matrix.
		FLA_Wilkshift_tridiag_opd( *d_last_m1,
		                           *e_last,
		                           *d_last,
    		                       &shift );

		// Perform a Francis step.
		r_val = FLA_Tevd_francis_v_opd_var1( m_A,
		                                     &shift,
		                                     g1,     rs_G,
		                                     buff_d, inc_d,
		                                     buff_e, inc_e );

		// Check for internal deflation.
		if ( r_val != FLA_SUCCESS )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_eigval_v_opt_var1: Internal deflation in col %d, eig %d\n", r_val, m_A - 1 );
			printf( "FLA_Tevd_eigval_v_opt_var1: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
			printf( "FLA_Tevd_eigval_v_opt_var1: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif

			// Set the off-diagonal element to zero.
			buff_e[ r_val*inc_e ] = 0.0;

			*n_iter = k + 1;
			return r_val;
		}

		/*------------------------------------------------------------*/
	}

	*n_iter = n_iter_allowed;
	return FLA_FAILURE;
}

