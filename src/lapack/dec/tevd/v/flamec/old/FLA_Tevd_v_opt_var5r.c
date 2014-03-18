
#include "FLAME.h"

FLA_Error FLA_Tevd_v_opt_var5r( FLA_Obj d, FLA_Obj e )
{
	FLA_Datatype datatype;
	int          m_A;
	int          inc_d;
	int          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	
	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			sgiv_t    rots;

			FLA_Tevd_v_ops_var5r( m_A,
			                      0,
			                      buff_d, inc_d,
			                      buff_e, inc_e,
			                      &rots );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dgiv_t    rots;

			FLA_Tevd_v_opd_var5r( m_A,
			                      0,
			                      buff_d, inc_d,
			                      buff_e, inc_e,
			                      &rots );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var5r( int       m_A,
                                int       ij,
                                float*    buff_d, int inc_d, 
                                float*    buff_e, int inc_e,
                                sgiv_t*   rots )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opd_var5r( int       m_A,
                                int       ij,
                                double*   buff_d, int inc_d, 
                                double*   buff_e, int inc_e,
                                dgiv_t*   rots )
{
	FLA_Error r_val;
	int       i;

	// Iterate from back to front until all that is left is a 2x2.
	for ( i = m_A - 1; i > 1; --i )
	{
		int m_ATL = i + 1;

		/*------------------------------------------------------------*/

		// Find an eigenvalue of ATL submatrix.
		r_val = FLA_Eigval_v_opd_var5( m_ATL,
		                               ij,
		                               buff_d, inc_d,
		                               buff_e, inc_e,
		                               rots );

		if ( r_val == FLA_FAILURE )
		{
			printf( "FLA_Tevd_v_opd_var5r: failed to converge (iteration %d)\n", i );
			FLA_Abort();
		}

#ifdef PRINTF
		printf( "FLA_Tevd_v_opd_var5r: found eig  %14.11f in col %3d (n=%d) after %3d iterations\n", buff_d[ r_val*inc_d ], r_val, m_ATL, -1 );
#endif

		// If r_val != i, then the eigenvalue converged elsewhere within A.
		// Therefore, we must recurse with subproblems.
		if ( r_val != i )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_v_opd_var5r: Internal deflation in col %d, eig %d\n", r_val, i );
			printf( "FLA_Tevd_v_opd_var5r: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
			printf( "FLA_Tevd_v_opd_var5r: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif

			//printf( "FLA_Tevd_v_opd_var5r: iter so far = %3d\n", k_total );

			// An eigenvalue converged somewhere within the diagonal (not at
			// either the end), so we have to recurse with two subproblems.
			{
				int       m_ATL = r_val;
				int       m_ABR = m_A - r_val - 1;
				int       ijT   = ij;
				int       ijB   = ij + m_ATL + 1;
				double*   dTL   = buff_d + (0      )*inc_d;
				double*   eTL   = buff_e + (0      )*inc_e;
				double*   dBR   = buff_d + (m_ATL+1)*inc_d;
				double*   eBR   = buff_e + (m_ATL+1)*inc_e;

				FLA_Tevd_v_opd_var5r( m_ATL,
				                      ijT,
				                      dTL, inc_d,
				                      eTL, inc_e,
				                      rots );
				FLA_Tevd_v_opd_var5r( m_ABR,
				                      ijB,
				                      dBR, inc_d,
				                      eBR, inc_e,
				                      rots );
				return FLA_SUCCESS;
			}
		}

		/*------------------------------------------------------------*/
	}

	//printf( "FLA_Tevd_v_opd_var5r: total iter  = %3d\n", k_total );

	// Skip 1x1 matrices (and submatrices) entirely.
	if ( m_A > 1 )
	{
		double* gamma1;
		double* sigma1;

/*
typedef struct fla_dgivens
{
	int        n_g;
	int        n_g_alloc;
	dcomplex** g;
	int*       j_off;
	int*       m_g;
} dgiv_t;
*/
		// Update the Givens structure.
		rots->n_g += 1;
		rots->g[ rots->n_g - 1 ]     = bl1_zallocv( 1 );
		rots->m_g[ rots->n_g - 1 ]   = 1;
		rots->j_off[ rots->n_g - 1 ] = ij;

		gamma1 = &(rots->g[ rots->n_g-1 ]->real);
		sigma1 = &(rots->g[ rots->n_g-1 ]->imag);

		// Find the eigenvalue decomposition of the remaining (or only) 2x2
		// submatrix.
		FLA_Hevv_2x2_opd( buff_d + (0  )*inc_d,
		                  buff_e + (0  )*inc_e,
		                  buff_d + (1  )*inc_d,
		                  buff_d + (0  )*inc_d,
		                  buff_d + (1  )*inc_d,
		                  gamma1,
		                  sigma1 );

#ifdef PRINTF
		printf( "FLA_Tevd_v_opd_var5r: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 1*inc_d ], 1, 2, 1 );
		printf( "FLA_Tevd_v_opd_var5r: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 0*inc_d ], 0, 2, 1 );
#endif
	}

	return FLA_SUCCESS;
}
