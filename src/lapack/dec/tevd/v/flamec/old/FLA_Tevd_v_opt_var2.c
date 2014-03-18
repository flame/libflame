
#include "FLAME.h"

FLA_Error FLA_Tevd_v_opt_var2( FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U )
{
	FLA_Datatype datatype;
	int          m_A, m_U, n_G;
	int          inc_d;
	int          inc_e;
	int          rs_G, cs_G;
	int          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A      = FLA_Obj_vector_dim( d );
	m_U      = FLA_Obj_length( U );
	n_G      = FLA_Obj_width( G );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	
	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_U     = FLA_Obj_row_stride( U );
	cs_U     = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_U = FLA_FLOAT_PTR( U );

			FLA_Tevd_v_ops_var2( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			FLA_Tevd_v_opd_var2( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			FLA_Tevd_v_opc_var2( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			FLA_Tevd_v_opz_var2( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_opd_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opc_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opz_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U )
{
	FLA_Error r_val;
	double    gamma, sigma;
	int       i, k;
	int       k_total = 0;
	int       k_weight = 0;

	// Iterate from back to front until all that is left is a 2x2.
	for ( i = m_A - 1; i > 1; --i )
	{
		int m_ATL = i + 1;

		/*------------------------------------------------------------*/

		// Find an eigenvalue of ATL submatrix.
		r_val = FLA_Tevd_eigval_v_opd_var1( m_ATL,
		                                    n_G,
		                                    buff_G, rs_G, cs_G,
		                                    buff_d, inc_d,
		                                    buff_e, inc_e,
		                                    &k );

		if ( r_val == FLA_FAILURE )
		{
			printf( "FLA_Tevd_v_opz_var2: failed to converge (iteration %d)\n", i );
			FLA_Abort();
		}

#ifdef PRINTF
		printf( "FLA_Tevd_v_opz_var2: found eig  %14.11f in col %3d (n=%d) after %3d iterations\n", buff_d[ r_val*inc_d ], r_val, m_ATL, k );
#endif

		// Apply the resulting Givens rotations to update the eigenvectors.
		FLA_Apply_G_rf_opz_var1( k,
		                         m_U,
		                         m_ATL,
		                         buff_G, rs_G, cs_G,
		                         buff_U, rs_U, cs_U );

		k_total += k;
		k_weight += i * k;

		// If r_val != i, then the eigenvalue converged elsewhere within A.
		// Therefore, we must recurse with subproblems.
		if ( r_val != i )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_v_opz_var2: Internal deflation in col %d, eig %d\n", r_val, i );
			printf( "FLA_Tevd_v_opz_var2: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
			printf( "FLA_Tevd_v_opz_var2: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif

			//printf( "FLA_Tevd_v_opz_var2: iter so far = %3d\n", k_total );

			// An eigenvalue converged somewhere within the diagonal (not at
			// either the end), so we have to recurse with two subproblems.
			{
				int       m_ATL = r_val;
				int       m_ABR = m_A - r_val - 1;
				double*   dTL   = buff_d + (0      )*inc_d;
				double*   eTL   = buff_e + (0      )*inc_e;
				dcomplex* GT    = buff_G + (0      )*rs_G;
				dcomplex* UL    = buff_U + (0      )*cs_U;
				double*   dBR   = buff_d + (m_ATL+1)*inc_d;
				double*   eBR   = buff_e + (m_ATL+1)*inc_e;
				dcomplex* GB    = buff_G + (m_ATL+1)*rs_G;
				dcomplex* UR    = buff_U + (m_ATL+1)*cs_U;

				FLA_Tevd_v_opz_var2( m_ATL,
				                     m_U,
				                     n_G,
				                     dTL, inc_d,
				                     eTL, inc_e,
				                     GT,  rs_G, cs_G,
				                     UL,  rs_U, cs_U );
				FLA_Tevd_v_opz_var2( m_ABR,
				                     m_U,
				                     n_G,
				                     dBR, inc_d,
				                     eBR, inc_e,
				                     GB,  rs_G, cs_G,
				                     UR,  rs_U, cs_U );
				return FLA_SUCCESS;
			}
		}

		/*------------------------------------------------------------*/
	}

	//printf( "FLA_Tevd_v_opz_var2: total iter  = %3d\n", k_total );

	// Skip 1x1 matrices (and submatrices) entirely.
	if ( m_A > 1 )
	{
		// Find the eigenvalue decomposition of the remaining (or only) 2x2
		// submatrix.
		FLA_Hevv_2x2_opd( buff_d + (0  )*inc_d,
		                  buff_e + (0  )*inc_e,
		                  buff_d + (1  )*inc_d,
		                  buff_d + (0  )*inc_d,
		                  buff_d + (1  )*inc_d,
		                  &gamma,
		                  &sigma );

		// Update the eigenvectors.
		FLA_Apply_G_mx2_opz( m_U,
		                     &gamma,
		                     &sigma,
		                     buff_U + (0  )*cs_U, rs_U,
		                     buff_U + (1  )*cs_U, rs_U );

#ifdef PRINTF
		printf( "FLA_Tevd_v_opz_var2: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 1*inc_d ], 1, 2, 1 );
		printf( "FLA_Tevd_v_opz_var2: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 0*inc_d ], 0, 2, 1 );
#endif
	}

	return FLA_SUCCESS;
}
