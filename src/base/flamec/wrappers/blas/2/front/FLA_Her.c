
#include "FLAME.h"

FLA_Error FLA_Her( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A )
{
	return FLA_Her_external( uplo, alpha, x, A );
}

