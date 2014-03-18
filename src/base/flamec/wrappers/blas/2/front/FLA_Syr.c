
#include "FLAME.h"

FLA_Error FLA_Syr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A )
{
	return FLA_Syr_external( uplo, alpha, x, A );
}

