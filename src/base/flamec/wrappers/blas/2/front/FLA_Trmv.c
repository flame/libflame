
#include "FLAME.h"

FLA_Error FLA_Trmv( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x )
{
	return FLA_Trmv_external( uplo, transa, diag, A, x );
}

