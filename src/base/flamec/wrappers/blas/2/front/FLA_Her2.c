
#include "FLAME.h"

FLA_Error FLA_Her2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
	return FLA_Her2_external( uplo, alpha, x, y, A );
}

