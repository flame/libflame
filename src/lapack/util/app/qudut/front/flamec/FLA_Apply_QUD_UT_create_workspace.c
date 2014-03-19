/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_create_workspace( FLA_Obj T, FLA_Obj R, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        m_W, n_W;

	datatype = FLA_Obj_datatype( T );
	m_W      = FLA_Obj_length( T );
	n_W      = FLA_Obj_width( R );

	FLA_Obj_create( datatype, m_W, n_W, 0, 0, W );

	return FLA_SUCCESS;
}

