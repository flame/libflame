/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm )
{
	FLA_Datatype datatype;
	integer          m_A;
	integer          inc_d;
	integer          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d    = FLA_FLOAT_PTR( d );
			float*    buff_e    = FLA_FLOAT_PTR( e );
			float*    buff_norm = FLA_FLOAT_PTR( norm );

			FLA_Norm1_tridiag_ops( m_A,
			                       buff_d, inc_d,
			                       buff_e, inc_e,
			                       buff_norm );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d    = FLA_DOUBLE_PTR( d );
			double*   buff_e    = FLA_DOUBLE_PTR( e );
			double*   buff_norm = FLA_DOUBLE_PTR( norm );

			FLA_Norm1_tridiag_opd( m_A,
			                       buff_d, inc_d,
			                       buff_e, inc_e,
			                       buff_norm );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Norm1_tridiag_ops( integer       m_A,
                                 float*    buff_d, integer inc_d, 
                                 float*    buff_e, integer inc_e,
                                 float*    norm )
{
	float*  d  = buff_d;
	float*  e  = buff_e;
	float   nm;
	integer     i;

	if ( m_A == 1 )
	{
		nm = fabs( *d );
	}
	else
	{
		float  d_first = d[ (0    )*inc_d ];
		float  e_first = e[ (0    )*inc_e ];
		float  e_last  = e[ (m_A-2)*inc_e ];
		float  d_last  = d[ (m_A-1)*inc_d ];

		// Record the maximum of the absolute row/column sums for the
		// first and last row/columns.
		nm = max( fabs( d_first ) + fabs( e_first ),
		          fabs( e_last  ) + fabs( d_last  ) );

		for ( i = 1; i < m_A - 2; ++i )
		{
			float  e0 = e[ (i-1)*inc_e ];
			float  e1 = e[ (i  )*inc_e ];
			float  d1 = d[ (i  )*inc_d ];

			// Update nm with the absolute row/column sum for the ith
			// row/column.
			nm = max( nm, fabs( e0 ) +
			              fabs( d1 ) +
			              fabs( e1 ) );
		}
	}

	*norm = nm;

	return FLA_SUCCESS;
}



FLA_Error FLA_Norm1_tridiag_opd( integer       m_A,
                                 double*   buff_d, integer inc_d, 
                                 double*   buff_e, integer inc_e,
                                 double*   norm )
{
	double* d  = buff_d;
	double* e  = buff_e;
	double  nm;
	integer     i;

	if ( m_A == 1 )
	{
		nm = fabs( *d );
	}
	else
	{
		double d_first = d[ (0    )*inc_d ];
		double e_first = e[ (0    )*inc_e ];
		double e_last  = e[ (m_A-2)*inc_e ];
		double d_last  = d[ (m_A-1)*inc_d ];

		// Record the maximum of the absolute row/column sums for the
		// first and last row/columns.
		nm = max( fabs( d_first ) + fabs( e_first ),
		          fabs( e_last  ) + fabs( d_last  ) );

		for ( i = 1; i < m_A - 2; ++i )
		{
			double e0 = e[ (i-1)*inc_e ];
			double e1 = e[ (i  )*inc_e ];
			double d1 = d[ (i  )*inc_d ];

			// Update nm with the absolute row/column sum for the ith
			// row/column.
			nm = max( nm, fabs( e0 ) +
			              fabs( d1 ) +
			              fabs( e1 ) );
		}
	}

	*norm = nm;

	return FLA_SUCCESS;
}

