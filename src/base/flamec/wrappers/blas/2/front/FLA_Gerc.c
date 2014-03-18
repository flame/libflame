
#include "FLAME.h"

FLA_Error FLA_Gerc( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
	return FLA_Gerc_external( conjx, conjy, alpha, x, y, A );
}

