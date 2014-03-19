/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_find_split( FLA_Obj d, FLA_Obj e )
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
			float*    buff_d      = FLA_FLOAT_PTR( d );
			float*    buff_e      = FLA_FLOAT_PTR( e );

			FLA_Bsvd_find_split_ops( m_A,
			                         buff_d, inc_d,
			                         buff_e, inc_e );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d      = FLA_DOUBLE_PTR( d );
			double*   buff_e      = FLA_DOUBLE_PTR( e );

			FLA_Bsvd_find_split_opd( m_A,
			                         buff_d, inc_d,
			                         buff_e, inc_e );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_split_ops( int       m_A,
                                   float*    buff_d, int inc_d, 
                                   float*    buff_e, int inc_e )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_split_opd( int       m_A,
                                   double*   buff_d, int inc_d, 
                                   double*   buff_e, int inc_e )
{
	int i;

	for ( i = 0; i < m_A - 1; ++i )
	{
		double* epsilon1 = buff_e + (i  )*inc_e;

		if ( *epsilon1 == 0.0 )
		{
			// Return index of split as i+1 since e_i is in the same
			// column as d_(i+1).
			return i + 1;
		}
	}

	// Return with no split found found.
	return FLA_FAILURE;
}

