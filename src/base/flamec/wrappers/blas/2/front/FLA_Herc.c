
#include "FLAME.h"

FLA_Error FLA_Herc( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A )
{
	return FLA_Herc_external( uplo, conj, alpha, x, A );
}

