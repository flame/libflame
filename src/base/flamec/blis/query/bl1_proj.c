
#include "blis1.h"

conj1_t bl1_proj_trans1_to_conj( trans1_t trans )
{
	if ( bl1_does_conj( trans ) ) return BLIS1_CONJUGATE;
	else                          return BLIS1_NO_CONJUGATE;
}

