
#include "FLAME.h"

FLA_Error FLA_Copyrt( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj B )
{
	return FLA_Copyrt_external( uplo, trans, A, B );
}

