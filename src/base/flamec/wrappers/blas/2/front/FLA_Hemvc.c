
#include "FLAME.h"

FLA_Error FLA_Hemvc( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
	return FLA_Hemvc_external( uplo, conja, alpha, A, x, beta, y );
}

