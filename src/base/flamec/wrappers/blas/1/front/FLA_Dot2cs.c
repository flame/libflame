
#include "FLAME.h"

FLA_Error FLA_Dot2cs( FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj beta, FLA_Obj rho )
{
	return FLA_Dot2cs_external( conj, alpha, x, y, beta, rho );
}

