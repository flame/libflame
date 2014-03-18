
#include "FLAME.h"

FLA_Error FLA_Trmvsx( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
	return FLA_Trmvsx_external( uplo, transa, diag, alpha, A, x, beta, y );
}

