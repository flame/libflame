/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_create_workspace( FLA_Obj T, FLA_Obj B, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        m_W, n_W;

	datatype = FLA_Obj_datatype( T );
	m_W      = FLA_Obj_length( T );
	n_W      = FLA_Obj_max_dim( B );

        FLA_Obj_create( datatype, m_W, n_W, 0, 0, W );

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_Q_UT_create_workspace_side( FLA_Side side, FLA_Obj T, FLA_Obj B, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        m_W, n_W;

	datatype = FLA_Obj_datatype( T );
	m_W      = FLA_Obj_length( T );

        if      ( side == FLA_LEFT  ) n_W = FLA_Obj_width( B );
        else if ( side == FLA_RIGHT ) n_W = FLA_Obj_length( B );
        else                          n_W = FLA_Obj_max_dim( B );

        FLA_Obj_create( datatype, m_W, n_W, 0, 0, W );

	return FLA_SUCCESS;
}

