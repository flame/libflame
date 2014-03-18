
#include "FLAME.h"

FLA_Error FLA_Ger( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
	return FLA_Ger_external( alpha, x, y, A );
}

