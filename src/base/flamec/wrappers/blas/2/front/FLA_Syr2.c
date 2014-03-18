
#include "FLAME.h"

FLA_Error FLA_Syr2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
	return FLA_Syr2_external( uplo, alpha, x, y, A );
}

