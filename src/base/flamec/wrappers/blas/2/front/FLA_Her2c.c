
#include "FLAME.h"

FLA_Error FLA_Her2c( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
	return FLA_Her2c_external( uplo, conj, alpha, x, y, A );
}

