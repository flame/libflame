
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

