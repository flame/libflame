
#include "FLAME.h"

FLA_Error FLA_Axpyrt( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
	return FLA_Axpyrt_external( uplo, trans, alpha, A, B );
}

