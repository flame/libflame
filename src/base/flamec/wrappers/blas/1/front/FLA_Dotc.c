
#include "FLAME.h"

FLA_Error FLA_Dotc( FLA_Conj conj, FLA_Obj x, FLA_Obj y, FLA_Obj rho )
{
	return FLA_Dotc_external( conj, x, y, rho );
}

