
#include "FLAME.h"

FLA_Error FLA_Gemvc( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
	return FLA_Gemvc_external( transa, conjx, alpha, A, x, beta, y );
}
