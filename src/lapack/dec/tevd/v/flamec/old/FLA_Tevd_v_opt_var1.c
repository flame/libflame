/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_v_opt_var1( FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U )
{
	FLA_Datatype datatype;
	integer          m_A, m_U, n_G;
	integer          inc_d;
	integer          inc_e;
	integer          rs_G, cs_G;
	integer          rs_U, cs_U;

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

			FLA_Tevd_v_ops_var1( m_A,
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

			FLA_Tevd_v_opd_var1( m_A,
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

			FLA_Tevd_v_opc_var1( m_A,
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

			FLA_Tevd_v_opz_var1( m_A,
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



FLA_Error FLA_Tevd_v_ops_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_U, integer rs_U, integer cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opd_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_U, integer rs_U, integer cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opc_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_U, integer rs_U, integer cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opz_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_U, integer rs_U, integer cs_U )
{
	double  gamma, sigma;
	integer     i, k;
	integer     k_total = 0;
	integer     k_weight = 0;

	for ( i = m_A - 1; i > 1; --i )
	{
		integer m_ATL = i + 1;

		/*------------------------------------------------------------*/

		// Find an eigenvalue of the top-left m_ATL-by-m_ATL matrix.
		FLA_Tevd_eigval_v_opd_var1( m_ATL,
		                            n_G,
		                            buff_G, rs_G, cs_G,
		                            buff_d, inc_d,
		                            buff_e, inc_e,
		                            &k );

		k_total += k;
		k_weight += i * k;
//printf( "FLA_Tevd_v_opz_var1: found eig  %14.11f in col %3d after %3d iterations\n", buff_d[ i*inc_d ], i, k );

		// Apply the Givens rotations to update the eigenvectors.
		FLA_Apply_G_rf_opz_var1( k,
		                         m_U,
		                         m_ATL,
		                         buff_G, rs_G, cs_G,
		                         buff_U, rs_U, cs_U );

		/*------------------------------------------------------------*/
	}


//printf( "FLA_Tevd_v_opz_var1: total iter:        %d\n", k_total );
//printf( "FLA_Tevd_v_opz_var1: weighted avg iter: %.3f\n", ( double ) k_weight / ( m_A * m_A / 2 ) );

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

	return FLA_SUCCESS;
}
