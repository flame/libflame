
#include "FLAME.h"

FLA_Error FLA_Axpys( FLA_Obj alpha0, FLA_Obj alpha1, FLA_Obj A, FLA_Obj beta, FLA_Obj B )
{
	return FLA_Axpys_external( alpha0, alpha1, A, beta, B );
}
